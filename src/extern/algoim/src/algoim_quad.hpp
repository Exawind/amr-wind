#ifndef ALGOIM_QUAD_H
#define ALGOIM_QUAD_H

// High-order accurate quadrature algorithms for implicitly defined domains in hyperrectangles, based on the
// algorithms developed in the paper:
//  - R. I. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in Hyperrectangles,
//   SIAM Journal on Scientific Computing, 37(2), A993-A1019 (2015), http://dx.doi.org/10.1137/140966290

#include <array>
#include <bitset>
#include <vector>
#include <algorithm>
#include "algoim_gaussquad.hpp"
#include "algoim_interval.hpp"
#include "algoim_boundingbox.hpp"
#include "algoim_multiloop.hpp"

namespace Algoim
{
    namespace detail
    {
        // Newton's method safeguarded by a standard bisection method. If the given function is monotone on the
        // interval [x0,x1] and the interval brackets a root, then there is only one root and the method is
        // guaranteed to find it (subject to the tolerance and the maximum number of steps).
        template<typename F>
        bool newtonBisection(const F& f, Real x0, Real x1, Real tol, int maxsteps, Real& root)
        {
            Real f0 = f(x0), f1 = f(x1);
            if ((f0 > 0.0 && f1 > 0.0) || (f0 < 0.0 && f1 < 0.0))
                return false;
            if (f0 == Real(0.0))
            {
                root = x0;
                return true;
            }
            if (f1 == Real(0.0))
            {
                root = x1;
                return true;
            }

            // x0 and x1 define the bracket; x0 always corresponds to negative value of f; x1 positive value of f
            if (f1 < 0.0)
                std::swap(x0, x1);

            // Initial guess is midpoint
            Real x = (x0 + x1)*0.5;
            Real fx, fpx;
            f.eval(x, fx, fpx);
            Real dx = x1 - x0;
            for (int step = 0; step < maxsteps; ++step)
            {
                if ((fpx*(x - x0) - fx)*(fpx*(x - x1) - fx) < 0.0 && fabs(fx) < fabs(dx*fpx)*0.5)
                {
                    // Step in Newton's method falls within bracket and is less than half the previous step size
                    dx = -fx / fpx;
                    Real xold = x;
                    x += dx;
                    if (xold == x)
                    {
                        root = x;
                        return true;
                    }
                }
                else
                {
                    // Revert to bisection
                    dx = (x1 - x0)*0.5;
                    x = x0 + dx;
                    if (x == x0)
                    {
                        root = x;
                        return true;
                    }
                }
                if (std::abs(dx) < tol)
                {
                    root = x;
                    return true;
                }
                f.eval(x, fx, fpx);
                if (fx == Real(0.0)) // Got very lucky
                {
                    root = x;
                    return true;
                }
                if (fx < 0.0)
                    x0 = x;
                else
                    x1 = x;
            }
            root = (x0 + x1)*0.5;
            return true;
        }

        // OneDimRootFind takes in a multi-dimensional functional, reimagines it as a one-dimensional function
        // by freezing all but one of the axes (with the values of the frozen coordinates are passed in the ctor),
        // and computes the roots of the resulting 1-D function in a given interval, assuming that the function
        // is monotone on that interval. Although a struct, OneDimRootFind is essentially a void function. 
        template<typename F, int N>
        struct OneDimRootFind
        {
            const F& f;
            mutable TinyVector<Real,N> x;
            const int dim;

            OneDimRootFind(const F& f, const TinyVector<Real,N>& x, int dim, Real xmin, Real xmax, std::vector<Real>& roots)
                : f(f), x(x), dim(dim)
            {
                Real t;
                if (newtonBisection(*this, xmin, xmax, 1e2*Real(std::numeric_limits<Real>::epsilon())*std::abs(xmax - xmin), 1024, t))
                    roots.push_back(t);
            }

            Real operator() (Real t) const
            {
                x(dim) = t;
                return f(x);
            }

            void eval(Real t, Real& fx, Real& fpx) const
            {
                x(dim) = t;
                fx = f(x);
                fpx = f.grad(x)(dim);
            }
        };

        // Computes the roots of a given multi-dimensional function when it is restricted to a given interval on
        // a given axis. isMonotone can be used to indicate whether the corresponding 1-D function is known to be
        // monotone, but can be set false if this is unknown. Note: this function internally uses Interval
        // arithemtic, and may invalidate prior values of Interval<N>::delta.
        template<typename F, int N>
        void rootFind(const F& f, const TinyVector<Real,N>& x, int dim, Real xmin, Real xmax, std::vector<Real>& roots, bool isMonotone, int level = 0)
        {
            // If the function is known (by the caller) to be monotone or the recursive level is sufficiently deep, call the monotone case.
            if (isMonotone || level >= 4)
            {
                OneDimRootFind<F,N>(f, x, dim, xmin, xmax, roots);
                return;
            }

            // The function is not known to be monotone on the given interval (although in a lot of practical cases it will be). Use
            // interval arithmetic to place bounds on the function's value and gradient.
            TinyVector<Interval<N>,N> xint;
            for (int i = 0; i < N; ++i)
            {
                if (i == dim)
                {
                    xint(i) = Interval<N>((xmin + xmax)*0.5, setComponent<Real,N>(Real(0.0), i, Real(1.0)));
                    Interval<N>::delta(i) = (xmax - xmin)*0.5;
                }
                else
                {
                    xint(i) = x(i);
                    Interval<N>::delta(i) = Real(0.0);
                }
            }
    
            // First test to see if the function's value is bounded away from zero; if it is, then there are no roots.
            Interval<N> val = f(xint);
            if (val.uniformSign())
                return;

            // Test to see if the function's derivative is bounded away from zero on the interval; if it is, then call the monotone case.
            TinyVector<Interval<N>,N> g = f.grad(xint);
            if (g(dim).uniformSign())
            {
                OneDimRootFind<F,N>(f, x, dim, xmin, xmax, roots);
                return;
            }

            // Tests were inconclusive -- divide the interval into two and try again.
            rootFind(f, x, dim, xmin, (xmin + xmax)*0.5, roots, false, level + 1);
            rootFind(f, x, dim, (xmin + xmax)*0.5, xmax, roots, false, level + 1);
        }

        // Determines the sign conditions for restricting a (possibly already restricted) level set function, i.e.,
        // sgn_L and sgn_U in [R. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and Volumes in
        // Hyperrectangles, SIAM J. Sci. Comput., Vol. 37, No. 2, pp. A993-A1019, http://dx.doi.org/10.1137/140966290].
        template<bool S>
        void determineSigns(bool positiveAbove, int sign, int& bottomSign, int& topSign)
        {
            if (S)
            {
                // When evaluating a surface integral, if the function is positive above the height function, then
                // the bottom side must be negative and the top side must be positive; if the function is positive
                // below the height function, then the bottom side must be positive and the top side must be negative.
                bottomSign = positiveAbove? -1 : 1;
                topSign = positiveAbove? 1 : -1;
            }
            else
            {
                // When evaluating a volume integral:
                //   If integrating over the positive part:
                //      if positive above the height function: bottom = +/-, top = +
                //      if positive below the height function: bottom = +, top = +/-
                //   If integrating over the negative part:
                //      if positive above the height function: bottom = -, top = +/-
                //      if positive below the height function: bottom = +/-, top = -
                //   If integrating over both the positive and negative part (i.e. unrestricted in sign), keep it alive
                if (sign == 1)
                {
                    bottomSign = positiveAbove? 0 : 1;
                    topSign = positiveAbove? 1 : 0;
                }
                else if (sign == -1)
                {
                    bottomSign = positiveAbove? -1 : 0;
                    topSign = positiveAbove? 0 : -1;
                }
                else
                    bottomSign = topSign = 0;
            }
        }
    } // namespace detail

    // PsiCode encodes sign information of restricted level set functions on particular sides of
    // a hyperrectangle in a packed array of bits. The first N bits encodes side information, the
    // N+1st bit is true iff the sign == 0, while the N+2nd bit stores the sign if sign != 0.
    template<int N>
    struct PsiCode
    {
        std::bitset<N + 2> bits;
        PsiCode() {}
        PsiCode(const TinyVector<int,N>& sides, int sign)
        {
            for (int dim = 0; dim < N; ++dim)
                bits[dim] = sides(dim) == 1;
            if (sign == 0)
                bits[N] = 1;
            else
                bits[N] = 0, bits[N+1] = sign == 1;
        }
        // Modify an existing code by restriction in a particular dimension.
        PsiCode(const PsiCode& i, int dim, int side, int sign) : bits(i.bits)
        {
            bits[dim] = side == 1;
            if (sign == 0)
                bits[N] = 1;
            else
                bits[N] = 0, bits[N+1] = sign == 1;
        }
        inline int side(int dim) const
        {
            return bits[dim] ? 1 : 0;
        }
        inline int sign() const
        {
            return bits[N] ? 0 : (bits[N+1] ? 1 : -1);
        }
    };

    // M-dimensional integral of an N-dimensional function restricted to given implicitly defined domains
    template<int M, int N, typename Phi, typename F, bool S>
    struct ImplicitIntegral
    {
        const Phi& phi;
        F& f;
        const TinyVector<bool,N> free;
        std::array<PsiCode<N>,1 << (N - 1)> psi;
        int psiCount;
        const BoundingBox<Real,N> xrange;
        const int p;
        int e0;
        TinyVector<Interval<N>,N> xint;

        // Prune the given set of functions by checking for the existence of the interface. If a function is
        // uniformly positive or negative and is consistent with specified sign, it can be removed. If a
        // function is uniformly positive or negative but inconsistent with specified sign, the domain of
        // integration is empty.
        bool prune()
        {
            for (int i = 0; i < psiCount; )
            {
                for (int dim = 0; dim < N; ++dim)
                    if (!free(dim))
                        xint(dim).alpha = xrange(psi[i].side(dim))(dim);
                Interval<N> res = phi(xint);
                if (res.uniformSign())
                {
                    if ((res.alpha >= 0.0 && psi[i].sign() >= 0) || (res.alpha <= 0.0 && psi[i].sign() <= 0))
                    {
                        --psiCount;
                        std::swap(psi[i], psi[psiCount]);
                    }
                    else
                        return false;
                }
                else
                    ++i;
            }
            return true;
        }

        // Gaussian quadrature for when the domain of integration is determined to be the entire M-dimensional cube.
        void tensorProductIntegral()
        {
            for (MultiLoop<M> i(0,p); i; ++i)
            {
                TinyVector<Real,N> x;
                Real w = 1.0;
                for (int dim = 0, k = 0; dim < N; ++dim)
                {
                    if (free(dim))
                    {
                        x(dim) = xrange.min(dim) + xrange.extent(dim) * GaussQuad::x(p, i(k));
                        w *= xrange.extent(dim) * GaussQuad::w(p, i(k));
                        ++k;
                    }
                }
                f.evalIntegrand(x, w);
            }
        }

        // Given x, valid in all free variables barring e0, root find in the e0 direction on each of the
        // implicit functions, and apply Gaussian quadrature to each segment. Weights are multiplied upon going
        // back up the tree of recursive calls.
        void evalIntegrand(TinyVector<Real,N> x, Real w) const
        {
            // Each thread has its own storage for computing roots. These are not corrupted by the recursive
            // chain of evalIntegrand() calls since each call is in a different templated namespace. Moreover,
            // if using OpenMP tasking, each task is assumed tied and so one thread can only execute "evalIntegrand"
            // from start to end uninterrupted.
            static thread_local std::vector<Real> roots;
            roots.clear();

            if (S)
            {
                // Surface integral
                assert(M == N && N >= 2); // x is thus valid in all variables except e0
                detail::rootFind<Phi,N>(phi, x, e0, xrange.min(e0), xrange.max(e0), roots, true);
                assert(roots.size() <= 1);
                for (int i = 0; i < static_cast<int>(roots.size()); ++i)
                {
                    x(e0) = roots[i];
                    TinyVector<Real,N> g = phi.grad(x);
                    f.evalIntegrand(x, mag(g)/std::abs(g(e0))*w);
                }
                return;
            }

            // Partition [xmin(e0), xmax(e0)] by roots found
            roots.push_back(xrange.min(e0));
            for (int i = 0; i < psiCount; ++i)
            {
                for (int dim = 0; dim < N; ++dim)
                    if (!free(dim))
                        x(dim) = xrange(psi[i].side(dim))(dim);
                // x is now valid in all variables except e0
                detail::rootFind<Phi,N>(phi, x, e0, xrange.min(e0), xrange.max(e0), roots, M > 1);
            }
            std::sort(roots.begin() + 1, roots.end());
            roots.push_back(xrange.max(e0));

            // In rare cases, degenerate segments can be found, filter out with a tolerance
            Real tol = 10.0*Real(std::numeric_limits<Real>::epsilon())*(xrange.max(e0) - xrange.min(e0));

            // Loop over segments of divided interval
            for (int i = 0; i < static_cast<int>(roots.size()) - 1; ++i)
            {
                if (roots[i+1] - roots[i] < tol)
                    continue;

                // Evaluate sign of phi within segment and check for consistency with psi
                bool okay = true;
                x(e0) = (roots[i] + roots[i+1])*0.5;
                for (int j = 0; j < psiCount && okay; ++j)
                {
                    for (int dim = 0; dim < N; ++dim)
                        if (!free(dim))
                            x(dim) = xrange(psi[j].side(dim))(dim);
                    okay &= phi(x) > 0.0 ? (psi[j].sign() >= 0) : (psi[j].sign() <= 0);
                }
                if (!okay)
                    continue;

                for (int j = 0; j < p; ++j)
                {
                    x(e0) = roots[i] + (roots[i+1] - roots[i]) * GaussQuad::x(p, j);
                    Real gw = (roots[i+1] - roots[i]) * GaussQuad::w(p, j);
                    f.evalIntegrand(x, w * gw);
                }
            }
        }

        // Main calling engine; parameters with underscores are copied upon entry but modified internally in the ctor
        ImplicitIntegral(const Phi& phi, F& f, const TinyVector<bool,N>& free, const std::array<PsiCode<N>,1 << (N-1)>& psi_, int psiCount_, const BoundingBox<Real,N>& xrange, int p, int level = 0)
            : phi(phi), f(f), free(free), psi(psi_), psiCount(psiCount_), xrange(xrange), p(p)
        {
            // For the one-dimensional base case, evaluate the bottom-level integral.
            if (M == 1)
            {
                for (int dim = 0; dim < N; ++dim)
                    if (free(dim))
                        e0 = dim;
                evalIntegrand(Real(0.0), Real(1.0));
                return;
            }

            // Establish interval bounds for prune() and remaining part of ctor.
            for (int dim = 0; dim < N; ++dim)
            {
                if (free(dim))
                {
                    xint(dim) = Interval<N>(xrange.midpoint(dim), setComponent<Real,N>(Real(0.0), dim, Real(1.0)));
                    Interval<N>::delta(dim) = xrange.extent(dim)*0.5;
                }
                else
                {
                    xint(dim) = Interval<N>(Real(0.0)); // xint(dim).delta will be set per psi function
                    Interval<N>::delta(dim) = Real(0.0);
                }
            }

            // Prune list of psi functions: if prune procedure returns false, then the domain of integration is empty.
            if (!prune())
                return;

            // If all psi functions were pruned, then the volumetric integral domain is the entire hyperrectangle.
            if (psiCount == 0)
            {
                if (!S)
                    tensorProductIntegral();
                return;
            }

            // Among all monotone height function directions, choose the one that makes the associated height function look as flat as possible.
            // This is a modification to the criterion presented in [R. Saye, High-Order Quadrature Methods for Implicitly Defined Surfaces and
            // Volumes in Hyperrectangles, SIAM J. Sci. Comput., Vol. 37, No. 2, pp. A993-A1019, http://dx.doi.org/10.1137/140966290].
            e0 = -1;
            Real max_quan = 0.0;
            for (int dim = 0; dim < N; ++dim)
                if (!free(dim))
                    xint(dim).alpha = xrange(psi[0].side(dim))(dim);
            TinyVector<Interval<N>,N> g = phi.grad(xint);
            for (int dim = 0; dim < N; ++dim)
            {
                if (free(dim) && std::abs(g(dim).alpha) > 1.001*g(dim).maxDeviation())
                {
                    Real quan = std::abs(g(dim).alpha) * xrange.extent(dim);
                    if (quan > max_quan)
                    {
                        max_quan = quan;
                        e0 = dim;
                    }
                }
            }

            // Check compatibility with all implicit functions whilst simultaneously constructing new implicit functions.
            std::array<PsiCode<N>,1 << (N-1)> newPsi;
            int newPsiCount = 0;
            for (int i = 0; i < psiCount; ++i)
            {
                // Evaluate gradient in an interval
                for (int dim = 0; dim < N; ++dim)
                    if (!free(dim))
                        xint(dim).alpha = xrange(psi[i].side(dim))(dim);
                TinyVector<Interval<N>,N> g = phi.grad(xint);

                // Determine if derivative in e0 direction is bounded away from zero.
                bool directionOkay = e0 != -1 && g(e0).uniformSign();
                if (!directionOkay)
                {
                    if (level < 16)
                    {
                        // Direction is not a good one, divide the domain into two along the biggest free extent
                        Real maxext = 0.0;
                        int ind = -1;
                        for (int dim = 0; dim < N; ++dim)
                        {
                            if (free(dim))
                            {
                                Real ext = xrange.extent(dim);
                                if (ext > maxext)
                                {
                                    maxext = ext;
                                    ind = dim;
                                }
                            }
                        }
                        assert(ind >= 0);
                        Real xmid = xrange.midpoint(ind);
                        ImplicitIntegral<M,N,Phi,F,S>(phi, f, free, psi, psiCount, BoundingBox<Real,N>(xrange.min(), setComponent(xrange.max(), ind, xmid)), p, level + 1);
                        ImplicitIntegral<M,N,Phi,F,S>(phi, f, free, psi, psiCount, BoundingBox<Real,N>(setComponent(xrange.min(), ind, xmid), xrange.max()), p, level + 1);
                        return;
                    }
                    else
                    {
                        // Halt subdivision because we have recursively subdivided too deep; evaluate level set functions at
                        // the centre of box and check compatibility with signs.
                        TinyVector<Real,N> xpoint = xrange.midpoint();
                        bool okay = true;
                        for (int j = 0; j < static_cast<int>(psi.size()) && okay; ++j)
                        {
                            for (int dim = 0; dim < N; ++dim)
                                if (!free(dim))
                                    xpoint(dim) = xrange(psi[j].side(dim))(dim);
                            okay &= phi(xpoint) >= 0.0 ? (psi[j].sign() >= 0) : (psi[j].sign() <= 0);
                        }
                        if (okay)
                        {
                            Real measure = 1.0;
                            for (int dim = 0; dim < N; ++dim)
                                if (free(dim))
                                    measure *= xrange.extent(dim);
                            if (S)
                                assert(M == N);
                            else
                                f.evalIntegrand(xpoint, measure);
                        }
                        return;
                    }
                }
                
                // Direction is okay - build restricted level set functions and determine the appropriate signs
                int bottomSign, topSign;
                detail::determineSigns<S>(g(e0).alpha > 0.0, psi[i].sign(), bottomSign, topSign);
                newPsi[newPsiCount++] = PsiCode<N>(psi[i], e0, 0, bottomSign);
                newPsi[newPsiCount++] = PsiCode<N>(psi[i], e0, 1, topSign);
                assert(newPsiCount <= 1 << (N - 1));
            }

            // Dimension reduction call
            assert(e0 != -1);
            ImplicitIntegral<M-1,N,Phi,ImplicitIntegral<M,N,Phi,F,S>,false>(phi, *this, setComponent(free, e0, false), newPsi, newPsiCount, xrange, p);
        }
    };

    // Partial specialisation on M=0 as a dummy base case for the compiler
    template<int N, typename Phi, typename F, bool S>
    struct ImplicitIntegral<0,N,Phi,F,S>
    {
        ImplicitIntegral(const Phi&, F&, const TinyVector<bool,N>&, const std::array<PsiCode<N>,1 << (N-1)>&, int, const BoundingBox<Real,N>&, int) {}
    };

    // QuadratureRule records quadrature points and weights as they are created by ImplicitIntegral
    template<int N>
    struct QuadratureRule
    {
        struct Node
        {
            TinyVector<Real,N> x;
            Real w;
            Node(const TinyVector<Real,N>& x, Real w) : x(x), w(w) {}
        };
        std::vector<Node> nodes;

        // evalIntegrand records quadrature nodes when it is invoked by ImplicitIntegral
        void evalIntegrand(const TinyVector<Real,N>& x, Real w)
        {
            nodes.push_back({x, w});
        }

        // Evaluate an integral applied to a given functional
        template<typename F, typename R = typename std::result_of<F(const TinyVector<Real,N>&)>::type>
        R operator()(const F& f) const
        {
            R sum = 0;
            for (const auto& pt : nodes)
                sum += f(pt.x) * pt.w;
            return sum;
        }

        // Sum of weights, i.e., the measure of the integrated domain
        Real sumWeights() const
        {
            Real sum = 0;
            for (const auto& pt : nodes)
                sum += pt.w;
            return sum;
        }
    };

    // Output a QuadratureScheme as an XML .vtp file for visualisation in ParaView or anything else that
    // supports XML VTK files
    template<typename T, int N>
    void outputQuadratureRuleAsVtpXML(const QuadratureRule<N>& q, T& stream)
    {
        static_assert(N == 2 || N == 3, "outputQuadratureRuleAsVtpXML only supports 2D and 3D quadrature schemes");
        stream << "<?xml version=\"1.0\"?>\n";
        stream << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        stream << "<PolyData>\n";
        stream << "<Piece NumberOfPoints=\"" << q.nodes.size() << "\" NumberOfVerts=\"" << q.nodes.size() << "\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\">\n";
        stream << "<Points>\n";
        stream << "  <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">";
        for (const auto& pt : q.nodes)
        {
            stream << pt.x(0) << ' ' << pt.x(1) << ' ';
            if (N == 3)
                stream << pt.x(2) << ' ';
        }
        stream << "</DataArray>\n";
        stream << "</Points>\n";
        stream << "<Verts>\n";
        stream << "  <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
        for (size_t i = 0; i < q.nodes.size(); ++i)
            stream << i << ' ';
        stream <<  "</DataArray>\n";
        stream << "  <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">";
        for (size_t i = 1; i <= q.nodes.size(); ++i)
            stream << i << ' ';
        stream << "</DataArray>\n";
        stream << "</Verts>\n";
        stream << "<PointData Scalars=\"w\">\n";
        stream << "  <DataArray type=\"Float32\" Name=\"w\" NumberOfComponents=\"1\" format=\"ascii\">";
        for (const auto& pt : q.nodes)
            stream << pt.w << ' ';
        stream << "</DataArray>\n";
        stream << "</PointData>\n";
        stream << "</Piece>\n";
        stream << "</PolyData>\n";
        stream << "</VTKFile>\n";
    };

    /* Main engine for generating high-order accurate quadrature schemes for implicitly defined domains in
       hyperrectangles. The level set function is given by the function object phi, the hyperrectangle is
       specified by xrange, and qo determines the number of quadrature points in the underlying 1D Gaussian
       quadrature schemes. Specifically,
       - phi is a user-defined function object which evaluates the level set function and its gradient. It
         must implement both
           template<typename T> operator() (const blitz::TinyVector<T,N>& x) const
         and
           template<typename T> grad(const blitz::TinyVector<T,N>& x) const
         In the simplest case, T = Real (= double) and the role of phi is to simply evaluate the level set
         function (e.g., return x(0)*x(0) + x(1)*x(1) - 1; for a unit circle) and its gradient (e.g., return
         blitz::TinyVector<double,2>(2.0*x(0), 2.0*x(1));). However, it is crucial that these two member
         functions be templated on T in order to enable a certain kind of interval arithmetic, similar to
         automatic differentiation. In essence, the interval arithmetic automatically computes a first-order
         Taylor series (with bounded remainder) of the given level set function, and uses that to make
         decisions concerning the existence of the interface and what direction to use when converting the
         implicitly defined geometry into the graph of an implicitly defined height function. This requirement
         on phi being able to correctly perform interval arithmetic places restrictions on the type of level
         set functions quadGen can be applied to.
       - xrange is a user-specified bounding box, indicating the extent of the hyperrectangle in N dimensions
         to which the quadrature algorithm is applied.
       - dim is used to specify the type of quadrature:
          - If dim < 0, compute a volumetric quadrature scheme, whose domain is implicitly defined by
            {phi < 0} intersected with xrange.
          - If dim == N, compute a curved surface quadrature scheme, whose domain is implicitly defined by
            {phi == 0} intersected with xrange.
          - If 0 <= dim && dim < N, compute a flat surface quadrature scheme for one of the sides of the
            hyperrectangle, i.e., {phi < 0}, intersected with xrange, intersected with the face
            {x(dim) == xrange(side)(dim)}.
       - side is used only when 0 <= dim && dim < N and specifies which side of the hyperrectangle to restrict
         to, either side == 0 or side == 1 for the “left” or “right” face, respectively (with normal pointing
         in the direction of the dim-th axis).
       - qo specifies the degree of the underlying one-dimensional Gaussian quadrature scheme and must satisfy
         1 <= qo && qo <= 10.
    */
    template<int N, typename F>
    QuadratureRule<N> quadGen(const F& phi, const BoundingBox<Real,N>& xrange, int dim, int side, int qo)
    {
        QuadratureRule<N> q;
        std::array<PsiCode<N>,1 << (N - 1)> psi;
        TinyVector<bool,N> free = true;
        if (0 <= dim && dim < N)
        {
            // Volume integral for one of the sides of a hyperrectangle (in dimensions N - 1)
            assert(side == 0 || side == 1);
            psi[0] = PsiCode<N>(setComponent<int,N>(0, dim, side), -1);
            free(dim) = false;
            ImplicitIntegral<N-1,N,F,QuadratureRule<N>,false>(phi, q, free, psi, 1, xrange, qo);
            // The quadrature method is given a restricted level set function to work with, but does not actually
            // initialise the dim-th component of each quadrature node's position. Do so now.
            for (auto& node : q.nodes)
                node.x(dim) = xrange(side)(dim);
        }
        else if (dim == N)
        {
            // Surface integral
            psi[0] = PsiCode<N>(0, -1);
            ImplicitIntegral<N,N,F,QuadratureRule<N>,true>(phi, q, free, psi, 1, xrange, qo);
        }
        else
        {
            // Volume integral in the full N dimensions
            psi[0] = PsiCode<N>(0, -1);
            ImplicitIntegral<N,N,F,QuadratureRule<N>,false>(phi, q, free, psi, 1, xrange, qo);
        }
        return q;
    }
} // namespace Algoim

#endif
