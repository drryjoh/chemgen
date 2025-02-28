```terminal
cd build_kokkos
cmake ../../../third_party/kokkos -DCMAKE_INSTALL_PREFIX=../install_kokkos -DKokkos_ENABLE_THREADS=ON -DKokkos_ENABLE_SERIAL=ON
make -j 8 install
```

```
g++ -o parallel_kokkos test_kokkos.cpp -I./install_kokkos/include -L./install_kokkos/lib -lkokkoscore -std=c++17 -pthread
```

```
sysctl -n hw.logicalcpu
```

```
export KOKKOS_NUM_THREADS=16  
```
