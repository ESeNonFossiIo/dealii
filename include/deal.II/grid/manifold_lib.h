// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__manifold_lib_h
#define dealii__manifold_lib_h


#include <deal.II/base/config.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/base/function.h>
#include <deal.II/base/function_parser.h>

DEAL_II_NAMESPACE_OPEN

/**
 * Manifold description for a spherical space coordinate system.
 *
 * You can use this Manifold object to describe any sphere, circle,
 * hypersphere or hyperdisc in two or three dimensions, both as a co-dimension
 * one manifold descriptor or as co-dimension zero manifold descriptor.
 *
 * The two template arguments match the meaning of the two template arguments
 * in Triangulation<dim, spacedim>, however this Manifold can be used to
 * describe both thin and thick objects, and the behavior is identical when
 * dim <= spacedim, i.e., the functionality of PolarManifold<2,3> is
 * identical to PolarManifold<3,3>.
 *
 * The two dimensional implementation of this class works by transforming
 * points to spherical coordinates, taking the average in that coordinate
 * system, and then transforming back the point to Cartesian coordinates. For
 * the three dimensional case, we use a simpler approach: we take the average
 * of the norm of the points, and use this value to shift the average point
 * along the radial direction. In order for this manifold to work correctly,
 * it cannot be attached to cells containing the center of the coordinate
 * system. This point is a singular point of the coordinate transformation,
 * and there taking averages does not make any sense.
 *
 * This class is used in step-1 and step-2 to describe the boundaries of
 * circles. Its use is also discussed in the results section of step-6.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim = dim>
class PolarManifold : public ChartManifold<dim, spacedim, spacedim>
{
public:
  /**
   * The Constructor takes the center of the spherical coordinates system.
   * This class uses the pull_back and push_forward mechanism to transform
   * from Cartesian to spherical coordinate systems, taking into account the
   * periodicity of base Manifold in two dimensions, while in three dimensions
   * it takes the middle point, and project it along the radius using the
   * average radius of the surrounding points.
   */
  PolarManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * Pull back the given point from the Euclidean space. Will return the polar
   * coordinates associated with the point @p space_point. Only used when
   * spacedim = 2.
   */
  virtual Point<spacedim>
  pull_back(const Point<spacedim> &space_point) const;

  /**
   * Given a point in the spherical coordinate system, this method returns the
   * Euclidean coordinates associated to the polar coordinates @p chart_point.
   * Only used when spacedim = 3.
   */
  virtual Point<spacedim>
  push_forward(const Point<spacedim> &chart_point) const;

  /**
   * Given a point in the spacedim dimensional Euclidean space, this
   * method returns the derivatives of the function $F$ that maps from
   * the polar coordinate system to the Euclidean coordinate
   * system. In other words, it is a matrix of size
   * $\text{spacedim}\times\text{spacedim}$.
   *
   * This function is used in the computations required by the
   * get_tangent_vector() function.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual
  DerivativeForm<1,spacedim,spacedim>
  push_forward_gradient(const Point<spacedim> &chart_point) const;

  /**
   * The center of the spherical coordinate system.
   */
  const Point<spacedim> center;
private:

  /**
   * Helper function which returns the periodicity associated with this
   * coordinate system, according to dim, chartdim, and spacedim.
   */
  static Tensor<1,spacedim> get_periodicity();
};


/**
 * Manifold description for a spherical space coordinate system.
 *
 * You can use this Manifold object to describe any sphere, circle,
 * hypersphere or hyperdisc in two or three dimensions, both as a co-dimension
 * one manifold descriptor or as co-dimension zero manifold descriptor.
 *
 * The two template arguments match the meaning of the two template arguments
 * in Triangulation<dim, spacedim>, however this Manifold can be used to
 * describe both thin and thick objects, and the behavior is identical when
 * dim <= spacedim, i.e., the functionality of SphericalManifold<2,3> is
 * identical to SphericalManifold<3,3>.
 *
 *  PolarManifold reflects the usual notion of polar coordinates but can be
 *  a bad choice in the case we are interested in the north or south pole.
 *  Consider for istance the couple of points \f$x_1=(1,\pi/3,0)\f$ and
 *  \f$x_2=(1,\pi/3,\pi)\f$.
 *  These two point are connented (using a PolarManifold) by means of the curve
 *  \$[
 *  \begin{align}
 *    s: [0,1]  & \rightarrow &  \mathbb S^3 \\
 *            t & \mapsto     &  (1,\pi/3,0) + (0,0,t\pi)
 *  \$]
 *  This curve is not a geodesic an probably it is not we would choose. A better
 *  courve would be the one passing through the North pole:
 *  \[
 *   s(t) = x_1 \cos(\alpha(t)) + \kappa \times x_1 \sin(\alpha(t)) +
 *   \kappa ( \kappa \cdot x_1) (1-\cos(\alpha(t))).
 *  \]
 *  where $\kappa = \frac{x_1 \times \x_2}{\Vert x_1 \times \x_2 \Vert}$
 *  and $\alpha(t) = t * \arccos(x_1 * x_2) $ for $t\in[0,1]$.
 *  Indeed, this is a geodesic and completely avoid the singularities in the poles.
 *
 *  @note This is a corollary of the so called \emph{Rodrigues' rotation formula}.
 *
 * @ingroup manifold
 *
 * @author Mauro Bardelloni, Luca Heltai, 2016
 */
template <int dim, int spacedim = dim>
class SphericalManifold : public Manifold<dim, spacedim>
{
public:
  /**
   * The Constructor takes the center of the spherical coordinates system.
   * This class uses the pull_back and push_forward mechanism to transform
   * from Cartesian to spherical coordinate systems, taking into account the
   * periodicity of base Manifold in two dimensions, while in three dimensions
   * it takes the middle point, and project it along the radius using the
   * average radius of the surrounding points.
   */
  SphericalManifold(const Point<spacedim> center = Point<spacedim>());

  /**
   * TODO:
   */
  virtual
  Point<spacedim>
  get_new_point(const Point<spacedim> &p1,
                const Point<spacedim> &p2,
                const double w) const;

  virtual
  Point<spacedim>
  get_new_point (const Quadrature<spacedim> &quad) const
  {
    Assert(quad.size() > 0,
           ExcMessage("Quadrature should have at least one point."));

    Assert(std::abs(std::accumulate(quad.get_weights().begin(), quad.get_weights().end(), 0.0)-1.0) < 1e-10,
           ExcMessage("The weights for the individual points should sum to 1!"));

    Point<spacedim> p = quad.point(0);
    double w = quad.weight(0);

    for (unsigned int i=1; i<quad.size(); ++i)
      {
        if ( w != 0 )
          p = get_new_point(p, quad.point(i) , w/(quad.weight(i) + w) );
        else
          p = quad.point(i);
        w += quad.weight(i);
      }

    return p;
  };

  /**
   * TODO:
   */
  virtual
  Tensor<1,spacedim>
  get_tangent_vector (const Point<spacedim> &x1,
                      const Point<spacedim> &x2) const;

  /**
   * TODO:
   */
  virtual
  Point<spacedim>
  project_to_manifold (const std::vector<Point<spacedim> > &vertices,
                       const Point<spacedim> &candidate) const;

  /**
   * The center of the spherical coordinate system.
   */
  const Point<spacedim> center;
};

/**
 * Cylindrical Manifold description.  In three dimensions, points are
 * transformed using a cylindrical coordinate system along the <tt>x-</tt>,
 * <tt>y-</tt> or <tt>z</tt>-axis (when using the first constructor of this
 * class), or an arbitrarily oriented cylinder described by the direction of
 * its axis and a point located on the axis.
 *
 * This class was developed to be used in conjunction with the @p cylinder or
 * @p cylinder_shell functions of GridGenerator. This function will throw an
 * exception whenever spacedim is not equal to three.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim = dim>
class CylindricalManifold : public Manifold<dim,spacedim>
{
public:
  /**
   * Constructor. Using default values for the constructor arguments yields a
   * cylinder along the x-axis (<tt>axis=0</tt>). Choose <tt>axis=1</tt> or
   * <tt>axis=2</tt> for a tube along the y- or z-axis, respectively. The
   * tolerance value is used to determine if a point is on the axis.
   */
  CylindricalManifold (const unsigned int axis = 0,
                       const double tolerance = 1e-10);

  /**
   * Constructor. If constructed with this constructor, the manifold described
   * is a cylinder with an axis that points in direction #direction and goes
   * through the given #point_on_axis. The direction may be arbitrarily
   * scaled, and the given point may be any point on the axis. The tolerance
   * value is used to determine if a point is on the axis.
   */
  CylindricalManifold (const Point<spacedim> &direction,
                       const Point<spacedim> &point_on_axis,
                       const double tolerance = 1e-10);

  /**
   * Compute new points on the CylindricalManifold. See the documentation of
   * the base class for a detailed description of what this function does.
   */
  virtual Point<spacedim>
  get_new_point(const Quadrature<spacedim> &quad) const;

protected:
  /**
   * The direction vector of the axis.
   */
  const Point<spacedim> direction;

  /**
   * An arbitrary point on the axis.
   */
  const Point<spacedim> point_on_axis;

private:
  /**
   * Helper FlatManifold to compute tentative midpoints.
   */
  FlatManifold<dim,spacedim> flat_manifold;

  /**
   * Relative tolerance to measure zero distances.
   */
  double tolerance;
};


/**
 * Manifold description derived from ChartManifold, based on explicit
 * Function<spacedim> and Function<chartdim> objects describing the
 * push_forward() and pull_back() functions.
 *
 * You can use this Manifold object to describe any arbitrary shape domain, as
 * long as you can express it in terms of an invertible map, for which you
 * provide both the forward expression, and the inverse expression.
 *
 * In debug mode, a check is performed to verify that the transformations are
 * actually one the inverse of the other.
 *
 * @ingroup manifold
 *
 * @author Luca Heltai, 2014
 */
template <int dim, int spacedim=dim, int chartdim=dim>
class FunctionManifold : public ChartManifold<dim, spacedim, chartdim>
{
public:
  /**
   * Explicit functions constructor. Takes a push_forward function of spacedim
   * components, and a pull_back function of @p chartdim components. See the
   * documentation of the base class ChartManifold for the meaning of the
   * optional @p periodicity argument.
   *
   * The tolerance argument is used in debug mode to actually check that the
   * two functions are one the inverse of the other.
   */
  FunctionManifold(const Function<chartdim> &push_forward_function,
                   const Function<spacedim> &pull_back_function,
                   const Tensor<1,chartdim> &periodicity=Tensor<1,chartdim>(),
                   const double tolerance=1e-10);

  /**
   * Expressions constructor. Takes the expressions of the push_forward
   * function of spacedim components, and of the pull_back function of @p
   * chartdim components. See the documentation of the base class
   * ChartManifold for the meaning of the optional @p periodicity argument.
   *
   * The strings should be the readable by the default constructor of the
   * FunctionParser classes. You can specify custom variable expressions with
   * the last two optional arguments. If you don't, the default names are
   * used, i.e., "x,y,z".
   *
   * The tolerance argument is used in debug mode to actually check that the
   * two functions are one the inverse of the other.
   */
  FunctionManifold(const std::string push_forward_expression,
                   const std::string pull_back_expression,
                   const Tensor<1,chartdim> &periodicity=Tensor<1,chartdim>(),
                   const typename FunctionParser<spacedim>::ConstMap = typename FunctionParser<spacedim>::ConstMap(),
                   const std::string chart_vars=FunctionParser<chartdim>::default_variable_names(),
                   const std::string space_vars=FunctionParser<spacedim>::default_variable_names(),
                   const double tolerance=1e-10,
                   const double h=1e-8);

  /**
   * If needed, we delete the pointers we own.
   */
  ~FunctionManifold();

  /**
   * Given a point in the @p chartdim coordinate system, uses the
   * push_forward_function to compute the push_forward of points in @p
   * chartdim space dimensions to @p spacedim space dimensions.
   */
  virtual Point<spacedim>
  push_forward(const Point<chartdim> &chart_point) const;

  /**
   * Given a point in the chartdim dimensional Euclidean space, this
   * method returns the derivatives of the function $F$ that maps from
   * the sub_manifold coordinate system to the Euclidean coordinate
   * system. In other words, it is a matrix of size
   * $\text{spacedim}\times\text{chartdim}$.
   *
   * This function is used in the computations required by the
   * get_tangent_vector() function. The default implementation calls
   * the get_gradient() method of the
   * FunctionManifold::push_forward_function() member class. If you
   * construct this object using the constructor that takes two string
   * expression, then the default implementation of this method uses a
   * finite difference scheme to compute the gradients(see the
   * AutoDerivativeFunction() class for details), and you can specify
   * the size of the spatial step size at construction time with the
   * @p h parameter.
   *
   * Refer to the general documentation of this class for more information.
   */
  virtual
  DerivativeForm<1,chartdim,spacedim>
  push_forward_gradient(const Point<chartdim> &chart_point) const;

  /**
   * Given a point in the spacedim coordinate system, uses the
   * pull_back_function to compute the pull_back of points in @p spacedim
   * space dimensions to @p chartdim space dimensions.
   */
  virtual Point<chartdim>
  pull_back(const Point<spacedim> &space_point) const;

private:
  /**
   * Constants for the FunctionParser classes.
   */
  const typename FunctionParser<spacedim>::ConstMap const_map;

  /**
   * Pointer to the push_forward function.
   */
  SmartPointer<const Function<chartdim>,
               FunctionManifold<dim,spacedim,chartdim> > push_forward_function;

  /**
   * Pointer to the pull_back function.
   */
  SmartPointer<const Function<spacedim>,
               FunctionManifold<dim,spacedim,chartdim> > pull_back_function;

  /**
   * Relative tolerance. In debug mode, we check that the two functions
   * provided at construction time are actually one the inverse of the other.
   * This value is used as relative tolerance in this check.
   */
  const double tolerance;

  /**
   * Check ownership of the smart pointers. Indicates whether this class is
   * the owner of the objects pointed to by the previous two member variables.
   * This value is set in the constructor of the class. If @p true, then the
   * destructor will delete the function objects pointed to be the two
   * pointers.
   */
  const bool owns_pointers;
};


/**
 * Manifold description for the surface of a Torus in three dimensions. The
 * Torus is assumed to be in the x-z plane. The reference coordinate system
 * is given by the angle $phi$ around the y axis, the angle $theta$ around
 * the centerline of the torus, and the distance to the centerline $w$
 * (between 0 and 1).
 *
 * This class was developed to be used in conjunction with
 * GridGenerator::torus.
 *
 * @ingroup manifold
 *
 * @author Timo Heister, 2016
 */
template <int dim>
class TorusManifold : public ChartManifold<dim,3,3>
{
public:
  static const int chartdim = 3;
  static const int spacedim = 3;

  /**
   * Constructor. Specify the radius of the centerline @p R and the radius
   * of the torus itself (@p r). The variables have the same meaning as
   * the parameters in GridGenerator::torus().
   */
  TorusManifold (const double R, const double r);

  /**
   * Pull back operation.
   */
  virtual Point<3>
  pull_back(const Point<3> &p) const;

  /**
   * Push forward operation.
   */
  virtual Point<3>
  push_forward(const Point<3> &chart_point) const;

  /**
   * Gradient.
   */
  virtual
  DerivativeForm<1,3,3>
  push_forward_gradient(const Point<3> &chart_point) const;

private:
  double r, R;
};

DEAL_II_NAMESPACE_CLOSE

#endif
