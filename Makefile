
CXX=clang++
# CXX=g++

LAPACK_L=-llapack
BLAS_L=-lblas
# OMP_L=-lomp
FLAGS=-O2 -g -std=c++14 -Wunused -Wall -Wshadow -Wno-sign-compare -fopenmp #-fopt-info #-DKAHAN
# FLAGS=-O3 -std=c++14 -Wunused -Wall -Wshadow -Wno-sign-compare -fopenmp #-DKAHAN

HEADERS=auxiliary.hh aux_math.hh grid_1d.hh integrate.hh polynomial.hh mo_1d.hh state_1d.hh system_1d.hh tensor.hh local_basis_1d.hh basis_1d.hh global_basis_1d.hh timing.hh system_3d.hh mo_3d.hh
OBJECTS=aux_math.o grid_1d.o integrate.o polynomial.o mo_1d.o state_1d.o system_1d.o tensor.o tensor-square.o basis_1d.o local_basis_1d.o global_basis_1d.o state_1d_4c2e.o state_1d_xci.o state_1d_apprx1.o
OBJECTS_P=$(patsubst %.o, build/%.o, $(OBJECTS))

build/%.o: %.cc $(HEADERS) Makefile
	$(CXX) $(FLAGS) -c $< -o $@

all: ho-1d # ho-3d

ho-1d: $(OBJECTS_P) Makefile ho-1d.cc
	$(CXX) $(FLAGS) ho-1d.cc $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L) # $(OMP_L)
# ho-3d: $(OBJECTS_P) Makefile ho-3d.cc
#	$(CXX) $(FLAGS) ho-3d.cc $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L) # $(OMP_L)

poisson-bench: Makefile poisson-bench.cc $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) poisson-bench.cc $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)



build/test-main.o: Makefile test-main.cc
	$(CXX) test-main.cc -c -o $@
test-misc: Makefile test-misc.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-misc.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-misc
# test-helmholtz: Makefile test-helmholtz.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
# 	$(CXX) $(FLAGS) test-helmholtz.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
# 	./test-helmholtz
# test-poisson: Makefile test-poisson.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
# 	$(CXX) $(FLAGS) test-poisson.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
# 	./test-poisson
# test-updateE-3d: Makefile test-updateE-3d.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
# 	$(CXX) $(FLAGS) test-updateE-3d.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
# 	./test-updateE-3d
test-integral-gaussian: Makefile test-integral-gaussian.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-integral-gaussian.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-integral-gaussian
test-integral-polynomial: Makefile test-integral-polynomial.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-integral-polynomial.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-integral-polynomial
test-tensor: Makefile test-tensor.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-tensor.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-tensor
test-dgemm: Makefile test-dgemm.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-dgemm.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-dgemm
test-grid: Makefile test-grid.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-grid.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-grid
test-scf-mo-norms: Makefile test-scf-mo-norms.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-scf-mo-norms.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-scf-mo-norms
test-xci-exact: Makefile test-xci-exact.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-xci-exact.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-xci-exact
test-cell-correlation-sym: Makefile test-cell-correlation-sym.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-cell-correlation-sym.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-cell-correlation-sym
test-cell-integrals: Makefile test-cell-integrals.cc build/test-main.o $(OBJECTS_P) $(HEADERS)
	$(CXX) $(FLAGS) test-cell-integrals.cc build/test-main.o $(OBJECTS_P) -o $@ $(BLAS_L) $(LAPACK_L)
	./test-cell-integrals

clean:
	rm -f test build/*o timings.log foo

distclean: clean
	rm -f ho-1d ho-3d test-integral-gaussian test-integral-polynomial test-tensor test-updateE-3d test-helmholtz poisson-bench

