#ifndef INCFLO_IF_LIST_
#define INCFLO_IF_LIST_

#include <AMReX_EB2.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <type_traits>

/********************************************************************************
 *                                                                              *
 * Union of a list (std::vector) of the same kind of implicit function.         *
 *                                                                              *
 ********************************************************************************/

template <class F>
class UnionListIF
{

public:
	UnionListIF(const amrex::Vector<F>& a_ifs)
		: m_ifs(a_ifs)
	{
		empty = !a_ifs.empty();
	}

	~UnionListIF()
	{
	}

	UnionListIF(const UnionListIF& rhs) = default;
	UnionListIF(UnionListIF&& rhs) = default;
	UnionListIF& operator=(const UnionListIF& rhs) = default;
	UnionListIF& operator=(UnionListIF&& rhs) = default;

	bool is_empty() const
	{
		return empty;
	}

	amrex::Real operator()(const amrex::RealArray& p) const
	{

		// NOTE: this assumes that m_ifs is not empty
		amrex::Real vmax = m_ifs[0](p);
		for(int i = 1; i < m_ifs.size(); i++)
		{
			amrex::Real vcur = m_ifs[i](p);
			if(vmax < vcur)
				vmax = vcur;
		}

		return vmax;

		// NOTE: this would have been nice, but for some reason it does not work :(
		// even though according to https://en.cppreference.com/w/cpp/algorithm/max
		// it should ... ?
		//F & f_max = std::max( m_ifs,
		//                      [&](const F & f1, const F & f2) {
		//                          return f1(p) < f2(p);
		//                      });
		//return f_max(p);
	}

private:
	amrex::Vector<F> m_ifs;
	bool empty;
};

/********************************************************************************
 *                                                                              *
 * Conditional Implicit Functions => CIF                                        *
 * Can be "turned off" based on run-time parameters                             *
 *                                                                              *
 ********************************************************************************/

/********************************************************************************
 *                                                                              *
 * Conditional implicit function: a "normal" implicit function is given the     *
 * additional property `active` (getter: `is_active()`, setter                  *
 * `set_active(bool)`) allowing ConditionalUnion and ConditionalIntersection    *
 * operations.                                                                  *
 *                                                                              *
 ********************************************************************************/

template <class F>
class CIF : public F
{

public:
	CIF(F&& f, bool a_active)
		: F(f)
		, m_active(a_active)
	{
	}

	~CIF()
	{
	}

	CIF(const CIF& rhs) = default;
	CIF(CIF&& rhs) = default;
	CIF& operator=(const CIF& rhs) = default;
	CIF& operator=(CIF&& rhs) = default;

	void set_active(bool a_active)
	{
		m_active = a_active;
	}
	bool is_active() const
	{
		return m_active;
	}

private:
	bool m_active;
};

/********************************************************************************
 *                                                                              *
 * Conditional Union of two (different) conditional implicit functions.         *
 * NOTE: at least one of the conditional implicit functions must be active.     *
 *                                                                              *
 ********************************************************************************/

template <class F1, class F2>
class UnionCIF
{

public:
	// Assuming that one always active
	UnionCIF(const F1& f1, const F2& f2)
		: m_f1(f1)
		, m_f2(f2)
		, f1_active(f1.is_acive())
		, f2_active(f2.is_active())
	{
		AMREX_ALWAYS_ASSERT_WITH_MESSAGE(f1.is_active() || f2.is_active(),
										 "You must specify at least one active implicit function.");
	}

	~UnionCIF()
	{
	}

	UnionCIF(const UnionCIF& rhs) = default;
	UnionCIF(UnionCIF&& rhs) = default;
	UnionCIF& operator=(const UnionCIF& rhs) = default;
	UnionCIF& operator=(UnionCIF&& rhs) = default;

	bool is_active() const
	{
		return true;
	}

	amrex::Real operator()(const amrex::RealArray& p) const
	{

		if(!f1_active)
			return m_f2(p);
		if(!f2_active)
			return m_f1(p);

		return std::max(m_f1(p), m_f2(p));
	}

private:
	F1 m_f1;
	F2 m_f2;
	bool f1_active;
	bool f2_active;
};

/********************************************************************************
 *                                                                              *
 * Conditional Intersection of two (different) conditional implicit functions.  *
 * NOTE: at least one of the conditional implicit functions must be active.     *
 *                                                                              *
 ********************************************************************************/

template <class F1, class F2>
class IntersectionCIF
{

public:
	// Assuming that one always active
	IntersectionCIF(const F1& f1, const F2& f2)
		: m_f1(f1)
		, m_f2(f2)
		, f1_active(f1.is_acive())
		, f2_active(f2.is_active())
	{
		AMREX_ALWAYS_ASSERT_WITH_MESSAGE(f1.is_active() || f2.is_active(),
										 "You must specify at least one active implicit function.");
	}

	~IntersectionCIF()
	{
	}

	IntersectionCIF(const IntersectionCIF& rhs) = default;
	IntersectionCIF(IntersectionCIF&& rhs) = default;
	IntersectionCIF& operator=(const IntersectionCIF& rhs) = default;
	IntersectionCIF& operator=(IntersectionCIF&& rhs) = default;

	bool is_active() const
	{
		return true;
	}

	amrex::Real operator()(const amrex::RealArray& p) const
	{

		if(!f1_active)
			return m_f2(p);
		if(!f2_active)
			return m_f1(p);

		return std::min(m_f1(p), m_f2(p));
	}

private:
	F1 m_f1;
	F2 m_f2;
	bool f1_active;
	bool f2_active;
};

#endif
