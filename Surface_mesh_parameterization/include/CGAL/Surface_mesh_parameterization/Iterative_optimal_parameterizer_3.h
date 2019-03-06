// Date          : 20 Feb 2019
// Author(s)     : Hardik Jain

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_OPTIMAL_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_OPTIMAL_PARAMETERIZER_3_H

#include <CGAL/license/Surface_mesh_parameterization.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Surface_mesh_parameterization/internal/angles.h>
#include <CGAL/Surface_mesh_parameterization/internal/kernel_traits.h>
#include <CGAL/Surface_mesh_parameterization/Error_code.h>

#include <CGAL/Surface_mesh_parameterization/Fixed_border_iterative_parameterizer_3.h>

#include <CGAL/Default.h>

#if defined(CGAL_EIGEN3_ENABLED)
#include <CGAL/Eigen_solver_traits.h>
#endif

#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/function_output_iterator.hpp>
#include <boost/unordered_set.hpp>
#include <CGAL/squared_distance_2.h> //for 2D functions
#include <CGAL/squared_distance_3.h> //for 3D functions
/// \file Iterative_optimal_parameterizer_3.h
#define LAMDA 1
namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMethods
///
/// The class `Iterative_optimal_parameterizer_3`
/// implements the *Optimal Parameterization* algorithm. This method
/// is sometimes called <i>OP</i> or just <i>Authalic parameterization</i> \cgalCite{cgal:dma-ipsm-02}.
///
/// A one-to-one mapping is guaranteed if the surface's border is mapped onto a convex polygon.
///
/// This class is a strategy  called by the main
/// parameterization algorithm `Fixed_border_iterative_parameterizer_3::parameterize()` and it:
/// - provides the template parameters `BorderParameterizer_` and `SolverTraits_`.
/// - implements `compute_w_ij()` to compute w_ij = (i, j), coefficient of the matrix A
///   for j neighbor vertex of i, based on Discrete Authalic Parameterization algorithm.
///
/// \cgalModels `Parameterizer_3`
///
/// \tparam TriangleMesh_ must be a model of `FaceGraph`.
///
/// \tparam BorderParameterizer_ is a Strategy to parameterize the surface border
///         and must be a model of `Parameterizer_3`.<br>
///         <b>%Default:</b>
/// \code
///   Circular_border_arc_length_parameterizer_3<TriangleMesh_>
/// \endcode
///
/// \tparam SolverTraits_ must be a model of `SparseLinearAlgebraTraits_d`.<br>
///         Note that the system is *not* symmetric because `Fixed_border_iterative_parameterizer_3`
///         does not remove border vertices from the system.<br>
///         <b>%Default:</b> If \ref thirdpartyEigen "Eigen" 3.1 (or greater) is available
///         and `CGAL_EIGEN3_ENABLED` is defined, then an overload of `Eigen_solver_traits`
///         is provided as default parameter:
/// \code
///   CGAL::Eigen_solver_traits<
///           Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
///                           Eigen::IncompleteLUT< double > > >
/// \endcode
///
/// \sa `CGAL::Surface_mesh_parameterization::Fixed_border_iterative_parameterizer_3<TriangleMesh, BorderParameterizer, SolverTraits>`
///
template < class TriangleMesh_,
class BorderParameterizer_ = Default,
class SolverTraits_ = Default>
class Iterative_optimal_parameterizer_3
    : public Fixed_border_iterative_parameterizer_3<
      TriangleMesh_,
      typename Default::Get<
      BorderParameterizer_,
      Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type,
      typename Default::Get<
      SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
      Eigen::IncompleteLUT<double> > > >::type >
#else
#pragma message("Error: You must either provide 'SolverTraits_' or link CGAL with the Eigen library")
SolverTraits_>::type > // no parameter provided, and Eigen is not enabled: don't compile
#endif
{
public:
#ifndef DOXYGEN_RUNNING
  typedef typename Default::Get<
      BorderParameterizer_,
      Circular_border_arc_length_parameterizer_3<TriangleMesh_> >::type  Border_parameterizer;

  typedef typename Default::Get<
      SolverTraits_,
#if defined(CGAL_EIGEN3_ENABLED)
      CGAL::Eigen_solver_traits<
      Eigen::BiCGSTAB<Eigen_sparse_matrix<double>::EigenType,
      Eigen::IncompleteLUT<double> > >
#else
  SolverTraits_ // no parameter provided, and Eigen is not enabled: so don't compile!
#endif
  >::type                                                     Solver_traits;
#else
  typedef Border_parameterizer_                               Border_parameterizer;
  typedef SolverTraits_                                       Solver_traits;
#endif

  typedef TriangleMesh_                                       TriangleMesh;

  // Private types
private:
  // Superclass
  typedef Fixed_border_iterative_parameterizer_3<TriangleMesh,
      Border_parameterizer,
      Solver_traits>         Base;

  // Private types
private:
  typedef typename boost::graph_traits<TriangleMesh>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;

  typedef CGAL::Vertex_around_target_circulator<TriangleMesh> vertex_around_target_circulator;
  typedef CGAL::Face_around_target_circulator<TriangleMesh> face_around_target_circulator;
  typedef CGAL::Vertex_around_face_circulator<TriangleMesh> vertex_around_face_circulator;

  // Traits subtypes:
  typedef typename Base::PPM                                   PPM;
  typedef typename Base::Kernel                                Kernel;
  typedef typename Base::NT                                    NT;
  typedef typename Base::Point_3                               Point_3;
  typedef typename Base::Point_2                               Point_2;
  typedef typename Base::Vector_3                              Vector_3;
  typedef boost::unordered_set<vertex_descriptor> Vertex_set;
  // Solver traits subtypes:
  typedef typename Solver_traits::Vector                       Vector;
  typedef typename Solver_traits::Matrix                       Matrix;

  typedef CGAL::dynamic_face_property_t<double>                                     Face_double_tag;
  typedef typename boost::property_map<TriangleMesh, Face_double_tag>::type         Face_Double_map;
  typedef CGAL::dynamic_vertex_property_t<double>                                   Vertex_double_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_double_tag>::const_type Vertex_Double_map;
  typedef CGAL::dynamic_vertex_property_t<int>                                      Vertex_int_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_int_tag>::const_type    Vertex_Int_map;
  typedef CGAL::dynamic_vertex_property_t<Point_2>                                  Vertex_point2_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_point2_tag>::const_type Vertex_point2_map;

  // Public operations
public:
  /// Constructor
  Iterative_optimal_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
      ///< %Object that maps the surface's border to 2D space.
      Solver_traits sparse_la = Solver_traits())
///< Traits object to access a sparse linear system.
: Fixed_border_iterative_parameterizer_3<TriangleMesh,
  Border_parameterizer,
  Solver_traits>(border_param, sparse_la),
  borderLength_3D(0)
  {

  }

  // Default copy constructor and operator =() are fine

  // Protected operations
protected:
  virtual NT compute_faceArea(TriangleMesh& mesh) {
    return 0;
  }

  virtual NT compute_borderLength_3D(TriangleMesh& mesh)  {
    borderLength_3D = Polygon_mesh_processing::longest_border(mesh).second;
    return 1.0;
  }

  virtual NT compute_faceWise_L2(TriangleMesh& mesh, Vertex_point2_map &uvmap)  {
    return 0;
  }

  virtual NT compute_vertexWise_L2(TriangleMesh& mesh, Vertex_set& vertices)  {
    return 0;
  }

  /// Compute the line i of matrix A for i inner vertex:
  virtual Error_code setup_inner_vertex_relations(Matrix& A,
      Matrix& A_prev,
      Vector&,
      Vector&,
      const TriangleMesh& mesh,
      vertex_descriptor vertex,
      Vertex_Int_map& vimap)
  {
    int i = get(vimap,vertex);

    // circulate over vertices around 'vertex' to compute w_ii and w_ijs
    NT w_ii = 0;
    int vertexIndex = 0;

    vertex_around_target_circulator v_j(halfedge(vertex, mesh), mesh), end = v_j;
    CGAL_For_all(v_j, end){
      // Call to virtual method to do the actual coefficient computation
      NT w_ij = -1.0 * compute_w_ij(mesh, vertex, v_j);
      if(w_ij > 0)
        w_ij *= -1.0;
      // w_ii = - sum of w_ijs
      w_ii -= w_ij;

      // Get j index
      int j = get(vimap, *v_j);

      // Set w_ij in matrix
      A.set_coef(i,j, w_ij, true /*new*/);
      A_prev.set_coef(i,j, w_ij, true);
      vertexIndex++;
    }

    if (vertexIndex < 2)
      return ERROR_NON_TRIANGULAR_MESH;

    // Set w_ii in matrix
    A.set_coef(i,i, w_ii, true /*new*/);
    return OK;
  }

  virtual double compute_sig_ij(TriangleMesh& mesh, Vertex_point2_map &uvmap, vertex_descriptor v_i, vertex_descriptor v_j, double& gamma) {
    //return (nedgeLength_2D(uvmap, v_i, v_j)/nedgeLength_3D(mesh, v_i, v_j));
    return pow(nedgeLength_3D(mesh, v_i, v_j)/nedgeLength_2D(uvmap, v_i, v_j),gamma);
  }

  virtual double distError(TriangleMesh& mesh)  {
    return 0;
  }

private:
  double borderLength_3D;
  /// Compute w_ij = (i, j), coefficient of matrix A for j neighbor vertex of i.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  NT compute_w_ij(const TriangleMesh& mesh,
      vertex_descriptor main_vertex_v_i,
      vertex_around_target_circulator neighbor_vertex_v_j) const
  {
    const PPM ppmap = get(vertex_point, mesh);

    const Point_3& position_v_i = get(ppmap, main_vertex_v_i);
    const Point_3& position_v_j = get(ppmap, *neighbor_vertex_v_j);

    // Compute the square norm of v_j -> v_i vector
    Vector_3 edge = position_v_i - position_v_j;
    double square_len = edge*edge;

    // Compute cotangent of (v_k,v_j,v_i) corner (i.e. cotan of v_j corner)
    // if v_k is the vertex before v_j when circulating around v_i
    vertex_around_target_circulator previous_vertex_v_k = neighbor_vertex_v_j;
    previous_vertex_v_k--;
    const Point_3& position_v_k = get(ppmap, *previous_vertex_v_k);
    NT cotg_psi_ij = internal::cotangent<Kernel>(position_v_k, position_v_j, position_v_i);
    NT cotg_beta_ij = internal::cotangent<Kernel>(position_v_i, position_v_k, position_v_j);

    // Compute cotangent of (v_i,v_j,v_l) corner (i.e. cotan of v_j corner)
    // if v_l is the vertex after v_j when circulating around v_i
    vertex_around_target_circulator next_vertex_v_l = neighbor_vertex_v_j;
    next_vertex_v_l++;
    const Point_3& position_v_l = get(ppmap,*next_vertex_v_l);
    NT cotg_theta_ij = internal::cotangent<Kernel>(position_v_i, position_v_j, position_v_l);
    NT cotg_alpha_ij = internal::cotangent<Kernel>(position_v_j, position_v_l, position_v_i);

    NT weight = 0.0;
    CGAL_assertion(square_len != 0.0); // two points are identical!
    if(square_len != 0.0) {
      weight = LAMDA * (cotg_beta_ij + cotg_alpha_ij);
      weight += (1-LAMDA) * ((cotg_psi_ij + cotg_theta_ij) / square_len);
    }
    return weight;
  }

  double nedgeLength_3D(TriangleMesh& mesh, vertex_descriptor v_i, vertex_descriptor v_j) {
    Point_3 p_i = mesh.point(v_i);
    Point_3 p_j = mesh.point(v_j);
    // return normalised distance
    return sqrt(CGAL::squared_distance(p_i,p_j))/borderLength_3D;
  }

  double nedgeLength_2D(Vertex_point2_map &uvmap, vertex_descriptor v_i, vertex_descriptor v_j) {
    Point_2 p_i = get(uvmap, v_i);
    Point_2 p_j = get(uvmap, v_j);
    // return normalised distance
    return sqrt(CGAL::squared_distance(p_i,p_j))/4.0;
  }

};


} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_AUTHALIC_PARAMETERIZER_3_H
