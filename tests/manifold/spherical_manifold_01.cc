// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check SphericalManifold for get_new_point and get_tangent_vector issues.

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/manifold_lib.h>


int
main()
{
  initlog();

  // Center and radius of the Ball
  Point<2> center(0.0, 0.0);
  double radius = center.norm();

  {
    const SphericalManifold<2,2> manifold(center);

    Point<2> P1(1.0, 0.0);
    Point<2> P2(0.0, 1.0);

    Point<2> Q = manifold.get_new_point(P1, P2, .5);

    deallog << "=================================" << std::endl;;
    deallog << manifold.get_new_point(P1, P2, .125) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .25) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .375) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .625) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .75) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .875) << std::endl;
    deallog << "=================================" << std::endl;
    deallog << manifold.get_new_point(P1, Q, .25) << std::endl;
    deallog << manifold.get_new_point(P1, Q, .5) << std::endl;
    deallog << manifold.get_new_point(P1, Q, .75) << std::endl;
    deallog << manifold.get_new_point(P1, P2,.5) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .25) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .5) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .75) << std::endl;
    deallog << "=================================" << std::endl;
  }

  {
    const SphericalManifold<1,2> manifold(center);

    Point<2> P1(1.0, 0.0);
    Point<2> P2(0.0, 1.0);

    Point<2> Q = manifold.get_new_point(P1, P2, .5);

    deallog << "=================================" << std::endl;;
    deallog << manifold.get_new_point(P1, P2, .125) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .25) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .375) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .5) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .625) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .75) << std::endl;
    deallog << manifold.get_new_point(P1, P2, .875) << std::endl;
    deallog << "=================================" << std::endl;
    deallog << manifold.get_new_point(P1, Q, .25) << std::endl;
    deallog << manifold.get_new_point(P1, Q, .5) << std::endl;
    deallog << manifold.get_new_point(P1, Q, .75) << std::endl;
    deallog << manifold.get_new_point(P1, P2,.5) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .25) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .5) << std::endl;
    deallog << manifold.get_new_point(Q, P2, .75) << std::endl;
    deallog << "=================================" << std::endl;
  }

  return 0;
}



