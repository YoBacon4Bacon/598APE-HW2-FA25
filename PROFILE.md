# Profiling
```
make bench_matmul.exe
sudo perf record -g --call-graph=dwarf ./bench_matmul.exe 0 16
sudo perf report --stdio --no-children -g graph,0.5,caller
```
