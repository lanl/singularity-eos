//------------------------------------------------------------------------------
// © 2021-2023. Triad National Security, LLC. All rights reserved.  This
// program was produced under U.S. Government contract 89233218CNA000001
// for Los Alamos National Laboratory (LANL), which is operated by Triad
// National Security, LLC for the U.S.  Department of Energy/National
// Nuclear Security Administration. All rights in the program are
// reserved by Triad National Security, LLC, and the U.S. Department of
// Energy/National Nuclear Security Administration. The Government is
// granted for itself and others acting on its behalf a nonexclusive,
// paid-up, irrevocable worldwide license in this material to reproduce,
// prepare derivative works, distribute copies to the public, perform
// publicly and display publicly, and to permit others to do so.
//------------------------------------------------------------------------------

// JMM: This is a dummy file to provide a source file for cmake to
// compile when built without fortran.
//
// The main reason to do this is that when we build without fortran
// bindings, the library is completely header only. However WITH
// fortran bindings, it is not.
//
// cmake supports a build mode for header only libraries by marking
// those libraries INTERFACE. Without this flag, a few things will go
// wrong:
// 1. CMake will be unable to infer a proper linker
//    so one must be specified by hand.
// 2. no dynamic library file (e.g., "*.a") will be generated,
//    but dependencies that link against it, such as another library,
//    or in our case, tests, will look for it. This will lead to
//    a failure of the downstream library at link time.
// 
// On the other hand, a library with implementation files
// CANNOT be marked INTERFACE, as otherwise source files will
// not be compiled.
//
// Unfortunately, switching the INTERFACE tag on and off in cmake
// is very cumbersome, as if a library is marked INTERFACE
// all its dependencies must ALSO be slurped in with an INTERFACE
// flag. This introduces significant branching in the cmake code
// and a lot of builer plate. For now, then, I simply include this
// empty source file.
//
// In the future, we could include centralized code here that we
// always want compiled, such as template instantiations or
// convenience functions.
