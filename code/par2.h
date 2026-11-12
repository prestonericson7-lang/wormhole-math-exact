/* par2.h -- the A9's second core for one function of the process at a time (10-29; the machine firmware bare v2.12's
 * par_post / par_join system calls). par2(fn, a0, a1): core 1 runs fn(a1) on its own stack while the caller runs
 * fn(a0); it returns once both have ended. Where the calls fail (Linux, an older firmware, core 1 busy or lost) both
 * halves run here one after the other: the same work, the same results. KQ_PAR2=0 keeps everything on this core.
 *
 * fn on core 1 must be pure computation and memory: no system call (one ends it, flagged), no malloc, and nothing the
 * caller's half writes. CORE 1'S OWN TLS (10-29): its thread pointer is a block of its own, made once from the program's
 * TLS image (PT_TLS: .tdata copied, .tbss zeroed) -- every __thread buffer (the attention kernels' vfb, ab, aTk; the RoPE
 * cache) is per core, as it is per thread under OpenMP. Before, core 1 ran on the caller's block: two cores, one buffer.
 * KQ_PAR2TLS=0: the caller's block, as before. */
#ifndef PAR2_H
#define PAR2_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/syscall.h>
#if defined(__arm__)
#include <link.h>
#endif

#if defined(__arm__)
/* fn returns here on core 1: the end of its function (svc 0xF0011) */
__attribute__((naked, noinline, used)) static void par2_tramp(void)
{
    __asm__ volatile("ldr r7, =0x000F0011\n\tsvc #0\n\t1: b 1b\n\t.ltorg");
}
#endif
static int par2_state = -1;                                 /* -1 untried, 0 off, 1 on */
static uint8_t *par2_stack;
static uint64_t par2_runs, par2_alone;
static uint32_t par2_tp1;                                   /* core 1's thread pointer (0: the caller's) */
static __thread uint32_t par2_probe = 0x5A17C0DEu;          /* par2_tls_make's layout check */
#define PAR2_STACK (512u << 10)
#if defined(__arm__)
typedef struct { const ElfW(Phdr) *ph; uintptr_t base; } par2_tls_t;
static int par2_tls_cb(struct dl_phdr_info *info, size_t size, void *data)
{
    (void)size;
    for (int i = 0; i < info->dlpi_phnum; i++)
        if (info->dlpi_phdr[i].p_type == PT_TLS) {
            ((par2_tls_t *)data)->ph = &info->dlpi_phdr[i];
            ((par2_tls_t *)data)->base = (uintptr_t)info->dlpi_addr;
            return 1;
        }
    return 0;
}
/* ARM's TLS is variant I: the thread pointer at the TCB (2 words), the program's block at tp + roundup(8, p_align),
 * glibc's thread descriptor just below tp (copied too: only read, if at all, by pure code) */
static void par2_tls_make(void)
{
    const char *e = getenv("KQ_PAR2TLS");
    if (e && !atoi(e)) return;
    par2_tls_t ti = { NULL, 0 };
    dl_iterate_phdr(par2_tls_cb, &ti);
    const ElfW(Phdr) *ph = ti.ph;
    if (!ph || ph->p_memsz == 0) return;
    const uint32_t pa = ph->p_align > 1u ? (uint32_t)ph->p_align : 1u;
    const uint32_t al = pa > 64u ? pa : 64u;                /* tp1's alignment: at least the segment's */
    const uint32_t off = (8u + pa - 1u) & ~(pa - 1u), pre = 2048u;   /* the block at tp + roundup(8, p_align): the linker's TPOFF */
    uint8_t *blk = (uint8_t *)malloc(pre + off + ph->p_memsz + 2u * al);
    if (!blk) return;
    uint8_t *tp1 = (uint8_t *)(((uintptr_t)blk + pre + al - 1u) & ~(uintptr_t)(al - 1u));
    uint32_t tp0;
    __asm__ volatile("mrc p15, 0, %0, c13, c0, 3" : "=r"(tp0));
    memcpy(tp1 - pre, (const uint8_t *)(uintptr_t)tp0 - pre, pre + 8u);             /* the descriptor and the TCB */
    memcpy(tp1 + off, (const void *)(ti.base + (uintptr_t)ph->p_vaddr), ph->p_filesz);   /* .tdata's image */
    memset(tp1 + off + ph->p_filesz, 0, ph->p_memsz - ph->p_filesz);                /* .tbss */
    /* the layout checked, not assumed: the probe (initialized 0x5A17C0DE, never written) must lie in the block at its
     * TPOFF and read back from the copy; otherwise core 1 keeps the caller's block */
    const uintptr_t po = (uintptr_t)&par2_probe - (uintptr_t)tp0;
    if (po < off || po + sizeof par2_probe > off + ph->p_memsz || *(const uint32_t *)(const void *)(tp1 + po) != 0x5A17C0DEu) {
        fprintf(stderr, "par2: the TLS layout is not the one assumed (probe at tp+%lu, block at tp+%u, %u bytes) -- core 1 shares the caller's\n",
                (unsigned long)po, (unsigned)off, (unsigned)ph->p_memsz);
        free(blk);
        return;
    }
    par2_tp1 = (uint32_t)(uintptr_t)tp1;
    fprintf(stderr, "par2: core 1 has its own TLS block (%u bytes, %u from the image; probe at tp+%lu)\n",
            (unsigned)ph->p_memsz, (unsigned)ph->p_filesz, (unsigned long)po);
}
#endif
static inline int par2_on(void)
{
    if (par2_state < 0) {
        par2_state = 0;
#if defined(__arm__)
        const char *e = getenv("KQ_PAR2");
        if (!(e && !atoi(e)) && (par2_stack = (uint8_t *)malloc(PAR2_STACK))) { par2_state = 1; par2_tls_make(); }
#endif
    }
    return par2_state;
}
static inline void par2(void (*fn)(void *), void *a0, void *a1)
{
#if defined(__arm__)
    if (par2_on()) {
        uint32_t tls;
        __asm__ volatile("mrc p15, 0, %0, c13, c0, 3" : "=r"(tls));
        if (par2_tp1) tls = par2_tp1;                       /* core 1's own block */
        const long r = syscall(0xF0010, (long)(uintptr_t)fn, (long)(uintptr_t)a1,
                               (long)(uintptr_t)(par2_stack + PAR2_STACK - 64u), (long)(uintptr_t)par2_tramp, (long)tls);
        if (r == 0) {
            fn(a0);
            const long j = syscall(0xF0012);
            if (j == 0) { par2_runs++; return; }
            fprintf(stderr, "par2: core 1's half failed (%ld) -- stopping\n", j);
            exit(70);                                       /* never a half-done result */
        }
        par2_state = 0;                                     /* no core 1 here: both halves on this one from now on */
        fprintf(stderr, "par2: no second core for the process here (%ld) -- one core\n", r);
    }
#endif
    par2_alone++;
    fn(a0);
    fn(a1);
}
#endif
