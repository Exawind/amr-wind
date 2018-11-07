/* blitz/config.h.  Generated from config.h.in by configure.  */
/* blitz/config.h.in.  Generated from configure.ac by autoheader.  */


/******************************************************************************
 * config.h           Compiler language support flags
 *
 * This file was generated automatically when running the configure script.
 * You should rerun configure each time you switch compilers, install new
 * standard libraries, or change compiler versions.
 *
 */



/* Macro for declaring aligned variables */
#define ALIGN_VARIABLE(vartype,varname,alignment) vartype varname;

/* Enable dimensions with > 2^31 elements (NOT IMPLEMENTED) */
/* #undef FULLY64BIT */

/* define if bool is a built-in type */
#define HAVE_BOOL /**/

/* define if the Boost library is available */
/* #undef HAVE_BOOST */

/* Define to 1 if you have the <boost/mpi.hpp> header file. */
/* #undef HAVE_BOOST_MPI_HPP */

/* define if the Boost::Serialization library is available */
/* #undef HAVE_BOOST_SERIALIZATION */

/* define if the compiler has <climits> header */
#define HAVE_CLIMITS /**/

/* define if the compiler has complex<T> */
#define HAVE_COMPLEX /**/

/* define if the compiler has standard complex<T> functions */
#define HAVE_COMPLEX_FCNS /**/

/* define if the compiler has complex math functions */
#define HAVE_COMPLEX_MATH1 /**/

/* define if the compiler has more complex math functions */
/* #undef HAVE_COMPLEX_MATH2 */

/* define if complex math functions are in namespace std */
#define HAVE_COMPLEX_MATH_IN_NAMESPACE_STD /**/

/* define if the compiler supports const_cast<> */
#define HAVE_CONST_CAST /**/

/* Define to 1 if you have the <cstring> header file. */
#define HAVE_CSTRING 1

/* define if the compiler supports default template parameters */
#define HAVE_DEFAULT_TEMPLATE_PARAMETERS /**/

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* define if the compiler supports dynamic_cast<> */
#define HAVE_DYNAMIC_CAST /**/

/* define if the compiler handle computations inside an enum */
#define HAVE_ENUM_COMPUTATIONS /**/

/* define if the compiler handles (int) casts in enum computations */
#define HAVE_ENUM_COMPUTATIONS_WITH_CAST /**/

/* define if the compiler supports exceptions */
#define HAVE_EXCEPTIONS /**/

/* define if the compiler supports the explicit keyword */
#define HAVE_EXPLICIT /**/

/* define if the compiler supports explicit template function qualification */
#define HAVE_EXPLICIT_TEMPLATE_FUNCTION_QUALIFICATION /**/

/* define if the compiler recognizes the full specialization syntax */
#define HAVE_FULL_SPECIALIZATION_SYNTAX /**/

/* define if the compiler supports function templates with non-type parameters
   */
#define HAVE_FUNCTION_NONTYPE_PARAMETERS /**/

/* define if the compiler supports IEEE math library */
#define HAVE_IEEE_MATH /**/

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `m' library (-lm). */
/* #undef HAVE_LIBM */

/* Define to 1 if you have the `papi' library (-lpapi). */
/* #undef HAVE_LIBPAPI */

/* define if the compiler supports member constants */
#define HAVE_MEMBER_CONSTANTS /**/

/* define if the compiler supports member templates */
#define HAVE_MEMBER_TEMPLATES /**/

/* define if the compiler supports member templates outside the class
   declaration */
#define HAVE_MEMBER_TEMPLATES_OUTSIDE_CLASS /**/

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* define if the compiler supports the mutable keyword */
#define HAVE_MUTABLE /**/

/* define if the compiler implements namespaces */
#define HAVE_NAMESPACES /**/

/* define if the compiler supports the Numerical C Extensions Group restrict
   keyword */
/* #undef HAVE_NCEG_RESTRICT */

/* define if the compiler supports the __restrict__ keyword */
#define HAVE_NCEG_RESTRICT_EGCS /**/

/* define if the compiler has numeric_limits<T> */
#define HAVE_NUMERIC_LIMITS /**/

/* define if the compiler accepts the old for scoping rules */
/* #undef HAVE_OLD_FOR_SCOPING */

/* define if the compiler supports partial ordering */
#define HAVE_PARTIAL_ORDERING /**/

/* define if the compiler supports partial specialization */
#define HAVE_PARTIAL_SPECIALIZATION /**/

/* define if the compiler supports reinterpret_cast<> */
#define HAVE_REINTERPRET_CAST /**/

/* define if the compiler supports Run-Time Type Identification */
#define HAVE_RTTI /**/

/* define if the compiler has getrusage() function */
#define HAVE_RUSAGE /**/

/* define if the compiler supports static_cast<> */
#define HAVE_STATIC_CAST /**/

/* define if the compiler supports ISO C++ standard library */
#define HAVE_STD /**/

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* define if the compiler supports Standard Template Library */
#define HAVE_STL /**/

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* define if the compiler supports System V math library */
/* #undef HAVE_SYSTEM_V_MATH */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <tbb/atomic.h> header file. */
/* #undef HAVE_TBB_ATOMIC_H */

/* define if the compiler supports basic templates */
#define HAVE_TEMPLATES /**/

/* define if the compiler supports templates as template arguments */
#define HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS /**/

/* define if the compiler supports use of the template keyword as a qualifier
   */
#define HAVE_TEMPLATE_KEYWORD_QUALIFIER /**/

/* define if the compiler supports template-qualified base class specifiers */
#define HAVE_TEMPLATE_QUALIFIED_BASE_CLASS /**/

/* define if the compiler supports template-qualified return types */
#define HAVE_TEMPLATE_QUALIFIED_RETURN_TYPE /**/

/* define if the compiler supports function matching with argument types which
   are template scope-qualified */
#define HAVE_TEMPLATE_SCOPED_ARGUMENT_MATCHING /**/

/* define if the compiler recognizes typename */
#define HAVE_TYPENAME /**/

/* define if the compiler supports the vector type promotion mechanism */
#define HAVE_TYPE_PROMOTION /**/

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* define if the compiler supports numeric traits promotions */
#define HAVE_USE_NUMTRAIT /**/

/* define if the compiler has valarray<T> */
#define HAVE_VALARRAY /**/

/* define if the compiler has isnan function in namespace std */
#define ISNAN_IN_NAMESPACE_STD /**/

/* Define to the sub-directory where libtool stores uninstalled libraries. */
#define LT_OBJDIR ".libs/"

/* define if the compiler has C math abs(integer types) in namespace std */
#define MATH_ABSINT_IN_NAMESPACE_STD /**/

/* define if the compiler has C math functions in namespace std */
#define MATH_FN_IN_NAMESPACE_STD /**/

/* Name of package */
#define PACKAGE "blitz"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "blitz-support@lists.sourceforge.net"

/* Define to the full name of this package. */
#define PACKAGE_NAME "blitz"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "blitz 0.10"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "blitz"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "0.10"

/* Pad array lengths to SIMD width. */
/* #undef PAD_ARRAYS */

/* Set SIMD instruction width in bytes */
#define SIMD_WIDTH 1

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Enable Blitz thread-safety features */
/* #undef THREADSAFE */

/* Use TBB atomic types */
/* #undef THREADSAFE_USE_TBB */

/* Specifies whether compiler alignment pragmas should be used */
/* #undef USE_ALIGNMENT_PRAGMAS */

/* Version number of package */
#define VERSION "0.10"

/* CXX */
#define _compiler_name "c++"

/* CXXFLAGS */
#define _compiler_options ""

/* date */
#define _config_date "Mon Nov  5 10:21:49 GMT 2018"

/* uname -a */
#define _os_name "Linux knut-Precision-M4800 4.15.0-38-generic #41-Ubuntu SMP Wed Oct 10 10:59:38 UTC 2018 x86_64 x86_64 x86_64 GNU/Linux"

/* target */
#define _platform "x86_64-pc-linux-gnu"
