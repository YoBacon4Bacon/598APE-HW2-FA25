# Passing by reference
For `./bench_matmul.exe 0 16`, we see 36.69% of the time is spent
on `__memmove_evex_unaligned_erms` (and 7.59% on `__memset_evex_unaligned_erms`).
This is some sort of SIMD memory movement that I think is caused by passing
and returning `Poly`s by value.

Before:
```
Matrix size: 16x16, Mode: 0 (ct*pt)
ref_time_sec=0.000005, enc_time_sec=5.259056, rel_err=0.000000
```

After:
```
Matrix size: 16x16, Mode: 0 (ct*pt)
ref_time_sec=0.000006, enc_time_sec=3.299121, rel_err=0.000000
```
