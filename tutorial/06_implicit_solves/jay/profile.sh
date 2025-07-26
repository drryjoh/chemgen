export PATH="/Users/jsampa/homebrew/opt/binutils/bin:$PATH"
export CPPFLAGS="-I$(brew --prefix gperftools)/include"
export LDFLAGS="-L$(brew --prefix gperftools)/lib"

g++ $CPPFLAGS -O2 -g -fno-omit-frame-pointer \
  -I/Users/jsampa/homebrew/opt/yaml-cpp/include \
  -L/Users/jsampa/homebrew/opt/yaml-cpp/lib -lyaml-cpp \
  -I$CONDA_PREFIX/include/eigen3 \
  bin/src/chemgen.cpp $LDFLAGS -lprofiler -o profile
export CPUPROFILE=chemgen.prof
./profile
pprof --text --cum ./profile chemgen.prof > profile.txt 2>/dev/null #non-interactive
# pprof ./profile chemgen.prof #interactive