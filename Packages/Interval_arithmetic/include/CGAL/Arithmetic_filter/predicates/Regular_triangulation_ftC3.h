// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : include/CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC3.h
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ======================================================================

// This file is automatically generated by the script
// examples/Interval_arithmetic/filtered_predicate_converter.

#ifndef CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H
#define CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H

CGAL_BEGIN_NAMESPACE

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &ry,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &sx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &sy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &sz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &swt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &twt)
{
  CGAL_assertion(Interval_nt_advanced::want_exceptions);
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  try
  {
    Oriented_side result = power_testC3(
		px.interval(),
		py.interval(),
		pz.interval(),
		pwt.interval(),
		qx.interval(),
		qy.interval(),
		qz.interval(),
		qwt.interval(),
		rx.interval(),
		ry.interval(),
		rz.interval(),
		rwt.interval(),
		sx.interval(),
		sy.interval(),
		sz.interval(),
		swt.interval(),
		tx.interval(),
		ty.interval(),
		tz.interval(),
		twt.interval());
    FPU_set_cw(backup);
    return result;
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
    FPU_set_cw(backup);
    return power_testC3(
		px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		rx.exact(),
		ry.exact(),
		rz.exact(),
		rwt.exact(),
		sx.exact(),
		sy.exact(),
		sz.exact(),
		swt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

struct Static_Filtered_power_testC3_20
{
  static double _bound;
  static double _epsilon_0;
  // static unsigned number_of_failures; // ?

  // Call this function from the outside to update the context.
  static void new_bound (double b, error = 0)
  {
    _bound = b;
    // recompute the epsilons: "just" call it over Static_filter_error.
    // That's the tricky part that might not work for everything.
    (void) update_epsilons(b,b,b,b,b,b,_epsilon_0);
    // TODO: We should verify that all epsilons have really been updated.
  }

  static Oriented_side update_epsilon(
	const Static_filter_error &px,
	const Static_filter_error &py,
	const Static_filter_error &pz,
	const Static_filter_error &pwt,
	const Static_filter_error &qx,
	const Static_filter_error &qy,
	const Static_filter_error &qz,
	const Static_filter_error &qwt,
	const Static_filter_error &rx,
	const Static_filter_error &ry,
	const Static_filter_error &rz,
	const Static_filter_error &rwt,
	const Static_filter_error &sx,
	const Static_filter_error &sy,
	const Static_filter_error &sz,
	const Static_filter_error &swt,
	const Static_filter_error &tx,
	const Static_filter_error &ty,
	const Static_filter_error &tz,
	const Static_filter_error &twt,
	double & epsilon_0)
  {
    typedef Static_filter_error FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      FT drx = rx - tx;
      FT dry = ry - ty;
      FT drz = rz - tz;
      FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
      FT dsx = sx - tx;
      FT dsy = sy - ty;
      FT dsz = sz - tz;
      FT dst = square(dsx) + square(dsy) + square(dsz) - swt + twt;
  
      return Oriented_side( - Static_Filtered_sign_of_determinant4x4_16::update_epsilon(dpx, dpy, dpz, dpt,
  						   dqx, dqy, dqz, dqt,
  						   drx, dry, drz, drt,
  						   dsx, dsy, dsz, dst,
  		epsilon_0));
  }

  static Oriented_side epsilon_variant(
	const Restricted_double &px,
	const Restricted_double &py,
	const Restricted_double &pz,
	const Restricted_double &pwt,
	const Restricted_double &qx,
	const Restricted_double &qy,
	const Restricted_double &qz,
	const Restricted_double &qwt,
	const Restricted_double &rx,
	const Restricted_double &ry,
	const Restricted_double &rz,
	const Restricted_double &rwt,
	const Restricted_double &sx,
	const Restricted_double &sy,
	const Restricted_double &sz,
	const Restricted_double &swt,
	const Restricted_double &tx,
	const Restricted_double &ty,
	const Restricted_double &tz,
	const Restricted_double &twt,
	const double & epsilon_0)
  {
    typedef Restricted_double FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      FT drx = rx - tx;
      FT dry = ry - ty;
      FT drz = rz - tz;
      FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
      FT dsx = sx - tx;
      FT dsy = sy - ty;
      FT dsz = sz - tz;
      FT dst = square(dsx) + square(dsy) + square(dsz) - swt + twt;
  
      return Oriented_side( - Static_Filtered_sign_of_determinant4x4_16::epsilon_variant(dpx, dpy, dpz, dpt,
  						   dqx, dqy, dqz, dqt,
  						   drx, dry, drz, drt,
  						   dsx, dsy, dsz, dst,
  		epsilon_0));
  }
};

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &ry,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &sx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &sy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &sz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &swt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &twt)
{
  bool re_adjusted = false;
  const double SAF_bound = Static_Filtered_power_testC3_20::_bound;

  // Check the bounds.  All arguments must be <= SAF_bound.
  if (	fabs(px.value()) > SAF_bound ||
	fabs(py.value()) > SAF_bound ||
	fabs(pz.value()) > SAF_bound ||
	fabs(pwt.value()) > SAF_bound ||
	fabs(qx.value()) > SAF_bound ||
	fabs(qy.value()) > SAF_bound ||
	fabs(qz.value()) > SAF_bound ||
	fabs(qwt.value()) > SAF_bound ||
	fabs(rx.value()) > SAF_bound ||
	fabs(ry.value()) > SAF_bound ||
	fabs(rz.value()) > SAF_bound ||
	fabs(rwt.value()) > SAF_bound ||
	fabs(sx.value()) > SAF_bound ||
	fabs(sy.value()) > SAF_bound ||
	fabs(sz.value()) > SAF_bound ||
	fabs(swt.value()) > SAF_bound ||
	fabs(tx.value()) > SAF_bound ||
	fabs(ty.value()) > SAF_bound ||
	fabs(tz.value()) > SAF_bound ||
	fabs(twt.value()) > SAF_bound)
  {
re_adjust:
    // Compute the new bound.
    double NEW_bound = std::max(0.0, fabs(px.value()));
    NEW_bound = std::max(NEW_bound, fabs(py.value()));
    NEW_bound = std::max(NEW_bound, fabs(pz.value()));
    NEW_bound = std::max(NEW_bound, fabs(pwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(qx.value()));
    NEW_bound = std::max(NEW_bound, fabs(qy.value()));
    NEW_bound = std::max(NEW_bound, fabs(qz.value()));
    NEW_bound = std::max(NEW_bound, fabs(qwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(rx.value()));
    NEW_bound = std::max(NEW_bound, fabs(ry.value()));
    NEW_bound = std::max(NEW_bound, fabs(rz.value()));
    NEW_bound = std::max(NEW_bound, fabs(rwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(sx.value()));
    NEW_bound = std::max(NEW_bound, fabs(sy.value()));
    NEW_bound = std::max(NEW_bound, fabs(sz.value()));
    NEW_bound = std::max(NEW_bound, fabs(swt.value()));
    NEW_bound = std::max(NEW_bound, fabs(tx.value()));
    NEW_bound = std::max(NEW_bound, fabs(ty.value()));
    NEW_bound = std::max(NEW_bound, fabs(tz.value()));
    NEW_bound = std::max(NEW_bound, fabs(twt.value()));
    // Re-adjust the context.
    Static_Filtered_power_testC3_20::new_bound(NEW_bound);
  }

  try
  {
    return Static_Filtered_power_testC3_20::epsilon_variant(px.to_double(),
		py.to_double(),
		pz.to_double(),
		pwt.to_double(),
		qx.to_double(),
		qy.to_double(),
		qz.to_double(),
		qwt.to_double(),
		rx.to_double(),
		ry.to_double(),
		rz.to_double(),
		rwt.to_double(),
		sx.to_double(),
		sy.to_double(),
		sz.to_double(),
		swt.to_double(),
		tx.to_double(),
		ty.to_double(),
		tz.to_double(),
		twt.to_double(),
		SAF_epsilon_0);
  }
  catch (Restricted_double::unsafe_comparison)
  {
    // It failed, we re-adjust once.
    if (!re_adjusted) {
      re_adjusted = true;
      goto re_adjust;
    }
    // This scheme definitely fails => exact computation (filtered_exact<> ?).
    return power_testC3(px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		rx.exact(),
		ry.exact(),
		rz.exact(),
		rwt.exact(),
		sx.exact(),
		sy.exact(),
		sz.exact(),
		swt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &ry,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &rwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &twt)
{
  CGAL_assertion(Interval_nt_advanced::want_exceptions);
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  try
  {
    Oriented_side result = power_testC3(
		px.interval(),
		py.interval(),
		pz.interval(),
		pwt.interval(),
		qx.interval(),
		qy.interval(),
		qz.interval(),
		qwt.interval(),
		rx.interval(),
		ry.interval(),
		rz.interval(),
		rwt.interval(),
		tx.interval(),
		ty.interval(),
		tz.interval(),
		twt.interval());
    FPU_set_cw(backup);
    return result;
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
    FPU_set_cw(backup);
    return power_testC3(
		px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		rx.exact(),
		ry.exact(),
		rz.exact(),
		rwt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

struct Static_Filtered_power_testC3_16
{
  static double _bound;
  static double _epsilon_0,_epsilon_1,_epsilon_2,_epsilon_3,_epsilon_4,_epsilon_5;
  // static unsigned number_of_failures; // ?

  // Call this function from the outside to update the context.
  static void new_bound (double b, error = 0)
  {
    _bound = b;
    // recompute the epsilons: "just" call it over Static_filter_error.
    // That's the tricky part that might not work for everything.
    (void) update_epsilons(b,b,b,b,b,b,_epsilon_0,_epsilon_1,_epsilon_2,_epsilon_3,_epsilon_4,_epsilon_5);
    // TODO: We should verify that all epsilons have really been updated.
  }

  static Oriented_side update_epsilon(
	const Static_filter_error &px,
	const Static_filter_error &py,
	const Static_filter_error &pz,
	const Static_filter_error &pwt,
	const Static_filter_error &qx,
	const Static_filter_error &qy,
	const Static_filter_error &qz,
	const Static_filter_error &qwt,
	const Static_filter_error &rx,
	const Static_filter_error &ry,
	const Static_filter_error &rz,
	const Static_filter_error &rwt,
	const Static_filter_error &tx,
	const Static_filter_error &ty,
	const Static_filter_error &tz,
	const Static_filter_error &twt,
	double & epsilon_0,
	double & epsilon_1,
	double & epsilon_2,
	double & epsilon_3,
	double & epsilon_4,
	double & epsilon_5)
  {
    typedef Static_filter_error FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      FT drx = rx - tx;
      FT dry = ry - ty;
      FT drz = rz - tz;
      FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
      Sign cmp;
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::update_epsilon(dpx, dpy, dpt,
  		                 dqx, dqy, dqt,
  				 drx, dry, drt,
  		epsilon_0);
      if (cmp != ZERO)
  	return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(px-rx, py-ry,
  		                                          qx-rx, qy-ry,
  		epsilon_1));
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::update_epsilon(dpx, dpz, dpt,
  		                 dqx, dqz, dqt,
  				 drx, drz, drt,
  		epsilon_2);
      if (cmp != ZERO)
  	return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(px-rx, pz-rz,
  		                                          qx-rx, qz-rz,
  		epsilon_3));
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::update_epsilon(dpy, dpz, dpt,
  		                 dqy, dqz, dqt,
  				 dry, drz, drt,
  		epsilon_4);
      return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(py-ry, pz-rz,
  		                                      qy-ry, qz-rz,
  		epsilon_5));
  }

  static Oriented_side epsilon_variant(
	const Restricted_double &px,
	const Restricted_double &py,
	const Restricted_double &pz,
	const Restricted_double &pwt,
	const Restricted_double &qx,
	const Restricted_double &qy,
	const Restricted_double &qz,
	const Restricted_double &qwt,
	const Restricted_double &rx,
	const Restricted_double &ry,
	const Restricted_double &rz,
	const Restricted_double &rwt,
	const Restricted_double &tx,
	const Restricted_double &ty,
	const Restricted_double &tz,
	const Restricted_double &twt,
	const double & epsilon_0,
	const double & epsilon_1,
	const double & epsilon_2,
	const double & epsilon_3,
	const double & epsilon_4,
	const double & epsilon_5)
  {
    typedef Restricted_double FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      FT drx = rx - tx;
      FT dry = ry - ty;
      FT drz = rz - tz;
      FT drt = square(drx) + square(dry) + square(drz) - rwt + twt;
      Sign cmp;
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::epsilon_variant(dpx, dpy, dpt,
  		                 dqx, dqy, dqt,
  				 drx, dry, drt,
  		epsilon_0);
      if (cmp != ZERO)
  	return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(px-rx, py-ry,
  		                                          qx-rx, qy-ry,
  		epsilon_1));
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::epsilon_variant(dpx, dpz, dpt,
  		                 dqx, dqz, dqt,
  				 drx, drz, drt,
  		epsilon_2);
      if (cmp != ZERO)
  	return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(px-rx, pz-rz,
  		                                          qx-rx, qz-rz,
  		epsilon_3));
  
      
      cmp = Static_Filtered_sign_of_determinant3x3_9::epsilon_variant(dpy, dpz, dpt,
  		                 dqy, dqz, dqt,
  				 dry, drz, drt,
  		epsilon_4);
      return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(py-ry, pz-rz,
  		                                      qy-ry, qz-rz,
  		epsilon_5));
  }
};

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &ry,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &rwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &twt)
{
  bool re_adjusted = false;
  const double SAF_bound = Static_Filtered_power_testC3_16::_bound;

  // Check the bounds.  All arguments must be <= SAF_bound.
  if (	fabs(px.value()) > SAF_bound ||
	fabs(py.value()) > SAF_bound ||
	fabs(pz.value()) > SAF_bound ||
	fabs(pwt.value()) > SAF_bound ||
	fabs(qx.value()) > SAF_bound ||
	fabs(qy.value()) > SAF_bound ||
	fabs(qz.value()) > SAF_bound ||
	fabs(qwt.value()) > SAF_bound ||
	fabs(rx.value()) > SAF_bound ||
	fabs(ry.value()) > SAF_bound ||
	fabs(rz.value()) > SAF_bound ||
	fabs(rwt.value()) > SAF_bound ||
	fabs(tx.value()) > SAF_bound ||
	fabs(ty.value()) > SAF_bound ||
	fabs(tz.value()) > SAF_bound ||
	fabs(twt.value()) > SAF_bound)
  {
re_adjust:
    // Compute the new bound.
    double NEW_bound = std::max(0.0, fabs(px.value()));
    NEW_bound = std::max(NEW_bound, fabs(py.value()));
    NEW_bound = std::max(NEW_bound, fabs(pz.value()));
    NEW_bound = std::max(NEW_bound, fabs(pwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(qx.value()));
    NEW_bound = std::max(NEW_bound, fabs(qy.value()));
    NEW_bound = std::max(NEW_bound, fabs(qz.value()));
    NEW_bound = std::max(NEW_bound, fabs(qwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(rx.value()));
    NEW_bound = std::max(NEW_bound, fabs(ry.value()));
    NEW_bound = std::max(NEW_bound, fabs(rz.value()));
    NEW_bound = std::max(NEW_bound, fabs(rwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(tx.value()));
    NEW_bound = std::max(NEW_bound, fabs(ty.value()));
    NEW_bound = std::max(NEW_bound, fabs(tz.value()));
    NEW_bound = std::max(NEW_bound, fabs(twt.value()));
    // Re-adjust the context.
    Static_Filtered_power_testC3_16::new_bound(NEW_bound);
  }

  try
  {
    return Static_Filtered_power_testC3_16::epsilon_variant(px.to_double(),
		py.to_double(),
		pz.to_double(),
		pwt.to_double(),
		qx.to_double(),
		qy.to_double(),
		qz.to_double(),
		qwt.to_double(),
		rx.to_double(),
		ry.to_double(),
		rz.to_double(),
		rwt.to_double(),
		tx.to_double(),
		ty.to_double(),
		tz.to_double(),
		twt.to_double(),
		SAF_epsilon_0,
		SAF_epsilon_1,
		SAF_epsilon_2,
		SAF_epsilon_3,
		SAF_epsilon_4,
		SAF_epsilon_5);
  }
  catch (Restricted_double::unsafe_comparison)
  {
    // It failed, we re-adjust once.
    if (!re_adjusted) {
      re_adjusted = true;
      goto re_adjust;
    }
    // This scheme definitely fails => exact computation (filtered_exact<> ?).
    return power_testC3(px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		rx.exact(),
		ry.exact(),
		rz.exact(),
		rwt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Dynamic, Protected, CGAL_IA_CACHE> &twt)
{
  CGAL_assertion(Interval_nt_advanced::want_exceptions);
  FPU_CW_t backup = FPU_get_and_set_cw(FPU_cw_up);
  try
  {
    Oriented_side result = power_testC3(
		px.interval(),
		py.interval(),
		pz.interval(),
		pwt.interval(),
		qx.interval(),
		qy.interval(),
		qz.interval(),
		qwt.interval(),
		tx.interval(),
		ty.interval(),
		tz.interval(),
		twt.interval());
    FPU_set_cw(backup);
    return result;
  } 
  catch (Interval_nt_advanced::unsafe_comparison)
  {
    FPU_set_cw(backup);
    return power_testC3(
		px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

struct Static_Filtered_power_testC3_12
{
  static double _bound;
  static double _epsilon_0,_epsilon_1,_epsilon_2,_epsilon_3,_epsilon_4,_epsilon_5;
  // static unsigned number_of_failures; // ?

  // Call this function from the outside to update the context.
  static void new_bound (double b, error = 0)
  {
    _bound = b;
    // recompute the epsilons: "just" call it over Static_filter_error.
    // That's the tricky part that might not work for everything.
    (void) update_epsilons(b,b,b,b,b,b,_epsilon_0,_epsilon_1,_epsilon_2,_epsilon_3,_epsilon_4,_epsilon_5);
    // TODO: We should verify that all epsilons have really been updated.
  }

  static Oriented_side update_epsilon(
	const Static_filter_error &px,
	const Static_filter_error &py,
	const Static_filter_error &pz,
	const Static_filter_error &pwt,
	const Static_filter_error &qx,
	const Static_filter_error &qy,
	const Static_filter_error &qz,
	const Static_filter_error &qwt,
	const Static_filter_error &tx,
	const Static_filter_error &ty,
	const Static_filter_error &tz,
	const Static_filter_error &twt,
	double & epsilon_0,
	double & epsilon_1,
	double & epsilon_2,
	double & epsilon_3,
	double & epsilon_4,
	double & epsilon_5)
  {
    typedef Static_filter_error FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      Comparison_result cmp;
  
      
      cmp = CGAL::Static_Filtered_compare_2::update_epsilon(px, qx,
  		epsilon_0);
      if (cmp != EQUAL)
          return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(dpx, dpt, dqx, dqt,
  		epsilon_1));
  
      
      cmp = CGAL::Static_Filtered_compare_2::update_epsilon(py, qy,
  		epsilon_2);
      if (cmp != EQUAL)
          return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(dpy, dpt, dqy, dqt,
  		epsilon_3));
  
      
      cmp = CGAL::Static_Filtered_compare_2::update_epsilon(pz, qz,
  		epsilon_4);
      return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::update_epsilon(dpz, dpt, dqz, dqt,
  		epsilon_5));
  }

  static Oriented_side epsilon_variant(
	const Restricted_double &px,
	const Restricted_double &py,
	const Restricted_double &pz,
	const Restricted_double &pwt,
	const Restricted_double &qx,
	const Restricted_double &qy,
	const Restricted_double &qz,
	const Restricted_double &qwt,
	const Restricted_double &tx,
	const Restricted_double &ty,
	const Restricted_double &tz,
	const Restricted_double &twt,
	const double & epsilon_0,
	const double & epsilon_1,
	const double & epsilon_2,
	const double & epsilon_3,
	const double & epsilon_4,
	const double & epsilon_5)
  {
    typedef Restricted_double FT;
  
      
      FT dpx = px - tx;
      FT dpy = py - ty;
      FT dpz = pz - tz;
      FT dpt = square(dpx) + square(dpy) + square(dpz) - pwt + twt;
      FT dqx = qx - tx;
      FT dqy = qy - ty;
      FT dqz = qz - tz;
      FT dqt = square(dqx) + square(dqy) + square(dqz) - qwt + twt;
      Comparison_result cmp;
  
      
      cmp = CGAL::Static_Filtered_compare_2::epsilon_variant(px, qx,
  		epsilon_0);
      if (cmp != EQUAL)
          return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(dpx, dpt, dqx, dqt,
  		epsilon_1));
  
      
      cmp = CGAL::Static_Filtered_compare_2::epsilon_variant(py, qy,
  		epsilon_2);
      if (cmp != EQUAL)
          return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(dpy, dpt, dqy, dqt,
  		epsilon_3));
  
      
      cmp = CGAL::Static_Filtered_compare_2::epsilon_variant(pz, qz,
  		epsilon_4);
      return Oriented_side(cmp * Static_Filtered_sign_of_determinant2x2_4::epsilon_variant(dpz, dpt, dqz, dqt,
  		epsilon_5));
  }
};

#ifndef CGAL_CFG_NO_EXPLICIT_TEMPLATE_FUNCTION_ARGUMENT_SPECIFICATION
template < class CGAL_IA_CT, class CGAL_IA_ET, class CGAL_IA_CACHE >
#endif
/*  */
Oriented_side
power_testC3(
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &px,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &py,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &pwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qy,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &qwt,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tx,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &ty,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &tz,
    const Filtered_exact <CGAL_IA_CT, CGAL_IA_ET, Static, Protected, CGAL_IA_CACHE> &twt)
{
  bool re_adjusted = false;
  const double SAF_bound = Static_Filtered_power_testC3_12::_bound;

  // Check the bounds.  All arguments must be <= SAF_bound.
  if (	fabs(px.value()) > SAF_bound ||
	fabs(py.value()) > SAF_bound ||
	fabs(pz.value()) > SAF_bound ||
	fabs(pwt.value()) > SAF_bound ||
	fabs(qx.value()) > SAF_bound ||
	fabs(qy.value()) > SAF_bound ||
	fabs(qz.value()) > SAF_bound ||
	fabs(qwt.value()) > SAF_bound ||
	fabs(tx.value()) > SAF_bound ||
	fabs(ty.value()) > SAF_bound ||
	fabs(tz.value()) > SAF_bound ||
	fabs(twt.value()) > SAF_bound)
  {
re_adjust:
    // Compute the new bound.
    double NEW_bound = std::max(0.0, fabs(px.value()));
    NEW_bound = std::max(NEW_bound, fabs(py.value()));
    NEW_bound = std::max(NEW_bound, fabs(pz.value()));
    NEW_bound = std::max(NEW_bound, fabs(pwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(qx.value()));
    NEW_bound = std::max(NEW_bound, fabs(qy.value()));
    NEW_bound = std::max(NEW_bound, fabs(qz.value()));
    NEW_bound = std::max(NEW_bound, fabs(qwt.value()));
    NEW_bound = std::max(NEW_bound, fabs(tx.value()));
    NEW_bound = std::max(NEW_bound, fabs(ty.value()));
    NEW_bound = std::max(NEW_bound, fabs(tz.value()));
    NEW_bound = std::max(NEW_bound, fabs(twt.value()));
    // Re-adjust the context.
    Static_Filtered_power_testC3_12::new_bound(NEW_bound);
  }

  try
  {
    return Static_Filtered_power_testC3_12::epsilon_variant(px.to_double(),
		py.to_double(),
		pz.to_double(),
		pwt.to_double(),
		qx.to_double(),
		qy.to_double(),
		qz.to_double(),
		qwt.to_double(),
		tx.to_double(),
		ty.to_double(),
		tz.to_double(),
		twt.to_double(),
		SAF_epsilon_0,
		SAF_epsilon_1,
		SAF_epsilon_2,
		SAF_epsilon_3,
		SAF_epsilon_4,
		SAF_epsilon_5);
  }
  catch (Restricted_double::unsafe_comparison)
  {
    // It failed, we re-adjust once.
    if (!re_adjusted) {
      re_adjusted = true;
      goto re_adjust;
    }
    // This scheme definitely fails => exact computation (filtered_exact<> ?).
    return power_testC3(px.exact(),
		py.exact(),
		pz.exact(),
		pwt.exact(),
		qx.exact(),
		qy.exact(),
		qz.exact(),
		qwt.exact(),
		tx.exact(),
		ty.exact(),
		tz.exact(),
		twt.exact());
  }
}

CGAL_END_NAMESPACE

#endif // CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H
