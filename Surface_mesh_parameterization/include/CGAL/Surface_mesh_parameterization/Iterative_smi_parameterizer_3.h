// Date          : 20 Feb 2019
// Author(s)     : Hardik Jain

#ifndef CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_SMI_PARAMETERIZER_3_H
#define CGAL_SURFACE_MESH_PARAMETERIZATION_ITERATIVE_SMI_PARAMETERIZER_3_H

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
/// \file Iterative_smi_parameterizer_3.h
#define LAMDA 1

namespace CGAL {

namespace Surface_mesh_parameterization {

/// \ingroup  PkgSurfaceMeshParameterizationMethods
///
/// The class `Iterative_smi_parameterizer_3`
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
class Iterative_smi_parameterizer_3
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
  typedef typename boost::property_map<TriangleMesh, Face_double_tag>::type   Face_Double_map;
  typedef CGAL::dynamic_vertex_property_t<double>                                   Vertex_double_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_double_tag>::const_type Vertex_Double_map;
  typedef CGAL::dynamic_vertex_property_t<Point_2>                                   Vertex_point2_tag;
  typedef typename boost::property_map<TriangleMesh, Vertex_point2_tag>::const_type Vertex_point2_map;

  //typedef typename TriangleMesh::Property_map<halfedge_descriptor, Point_2> SM_UV_pmap;

  // Public operations
public:
  /// Constructor
  Iterative_smi_parameterizer_3(Border_parameterizer border_param = Border_parameterizer(),
      ///< %Object that maps the surface's border to 2D space.
      Solver_traits sparse_la = Solver_traits())
///< Traits object to access a sparse linear system.
: Fixed_border_iterative_parameterizer_3<TriangleMesh,
  Border_parameterizer,
  Solver_traits>(border_param, sparse_la)
//  ,areaMap(get(Face_double_tag(), mesh))
  {
  }
  // Default copy constructor and operator =() are fine

private:
  TriangleMesh mesh;
  Vertex_Double_map vL2Map;
  Face_Double_map areaMap;
  Face_Double_map fL2Map;
  // Protected operations
protected:

  virtual NT compute_faceWise_L2(TriangleMesh& mesh, Vertex_point2_map &uvmap) {
    fL2Map = get(Face_double_tag(), mesh);
    BOOST_FOREACH(face_descriptor fd, faces(mesh))  {
      Point_2 uv_points[3];
      Point_3 mesh_points[3];
      int i = 0;
      BOOST_FOREACH(vertex_descriptor vd, CGAL::vertices_around_face(halfedge(fd,mesh),mesh)) {
        uv_points[i] = get(uvmap, vd);
        mesh_points[i] = mesh.point(vd);
        i++;
      }
      put(fL2Map, fd, getfL2(mesh_points, uv_points));
    }

    return 1.0;
  }

  virtual NT compute_faceArea(TriangleMesh& mesh) {
    // compute face area for mesh
    //typename TriangleMesh::template Property_map<face_descriptor, double> areaMap;
    //areaMap = mesh.add_property_map<face_descriptor, double>("f:area").first;
    areaMap = get(Face_double_tag(), mesh);
    BOOST_FOREACH(face_descriptor fd, faces(mesh))
        put(areaMap, fd, Polygon_mesh_processing::face_area(fd,mesh));
    return 1.0;
  }

  /// Compute w_ij = (i, j), coefficient of matrix A for j neighbor vertex of i.
  ///
  /// \param mesh a triangulated surface.
  /// \param main_vertex_v_i the vertex of `mesh` with index `i`
  /// \param neighbor_vertex_v_j the vertex of `mesh` with index `j`
  virtual NT compute_w_ij(const TriangleMesh& mesh,
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


  // compute_faceWise_L2
  virtual double getfL2(Point_3 mesh_points[], Point_2 uv_points[]) const {
    const double A = getA(uv_points);
    Point_3 Ss = getSs(mesh_points, uv_points, A);
    Point_3 St = getSt(mesh_points, uv_points, A);

    double a = innerProduct(Ss,Ss);
    double c = innerProduct(St,St);
    return sqrt((a+c)/2.0);
  }

  double getA(Point_2 uv_points[]) const {
    double A = (((uv_points[1].x()-uv_points[0].x())*(uv_points[2].y()-uv_points[0].y()))-((uv_points[2].x()-uv_points[0].x())*(uv_points[1].y()-uv_points[0].y())))/2;
    if (A == 0.0)
      return 1.0;
    return A;
  }

  Point_3 getSs(Point_3 mesh_points[3], Point_2 uv_points[3],const double &A) const {
    double dt0 = uv_points[1].y()-uv_points[2].y();
    double dt1 = uv_points[2].y()-uv_points[0].y();
    double dt2 = uv_points[0].y()-uv_points[1].y();
    Point_3 Ss (
        (mesh_points[0].x()*dt0 + mesh_points[1].x()*dt1 + mesh_points[2].x()*dt2 )/(2.0*A),
        (mesh_points[0].y()*dt0 + mesh_points[1].y()*dt1 + mesh_points[2].y()*dt2 )/(2.0*A),
        (mesh_points[0].z()*dt0 + mesh_points[1].z()*dt1 + mesh_points[2].z()*dt2 )/(2.0*A));
    return Ss;
  }

  Point_3 getSt(Point_3 mesh_points[3], Point_2 uv_points[3],const double &A) const  {
    double ds0 = uv_points[2].x()-uv_points[1].x();
    double ds1 = uv_points[0].x()-uv_points[2].x();
    double ds2 = uv_points[1].x()-uv_points[0].x();
    Point_3 St (
        (mesh_points[0].x()*ds0 + mesh_points[1].x()*ds1 +mesh_points[2].x()*ds2)/(2.0*A),
        (mesh_points[0].y()*ds0 + mesh_points[1].y()*ds1 +mesh_points[2].y()*ds2)/(2.0*A),
        (mesh_points[0].z()*ds0 + mesh_points[1].z()*ds1 +mesh_points[2].z()*ds2)/(2.0*A));
    return St;
  }

  double innerProduct(const Point_3& pointA, const Point_3& pointB) const {
    return ((pointA.x())*(pointB.x()) + (pointA.y())*(pointB.y()) + (pointA.z())*(pointB.z()));
  }



  // for compute_vertexWise_L2
  virtual NT compute_vertexWise_L2(TriangleMesh& mesh,
      Vertex_set& vertices) {
    // update weights for vertices
    vL2Map = get(Vertex_double_tag(), mesh);
    BOOST_FOREACH(vertex_descriptor v, vertices)
        put(vL2Map, v, getvL2(mesh, v));
    return 1.0;
  }

  // compute sigma to improve weights over iterations
  double getvL2(const TriangleMesh& mesh, vertex_descriptor &vertex)  const {
    bool debug = false;
    if(vertex == 711111)
      debug = true;
    if(debug)
      std::cout << "\t" << vertex << "\t" << std::flush;
    halfedge_descriptor hf = halfedge(vertex, mesh);
    if(debug)
      std::cout << hf << "\t" << std::flush;
    face_around_target_circulator f_j(hf, mesh), end_f_j = f_j;
    double varphi = 0.0;
    double localArea = 0.0;
    int i=0;
    CGAL_For_all(f_j, end_f_j)  {
      if(debug)
        std::cout << *f_j << "\t" << std::flush;
      if(*f_j > mesh.number_of_faces())
        continue;
      varphi += get(fL2Map,*f_j)*get(areaMap,*f_j);
      localArea += get(areaMap,*f_j);
      i++;
    }

    if(debug) {
      std::cout << varphi << "\t" << std::flush;
      std::cout << localArea << "\t" << std::flush;
    }

    if(mesh.is_border(vertex) && mesh.degree(vertex) != i+1)
      std::cerr << std::endl;
    else if(!mesh.is_border(vertex) && mesh.degree(vertex) != i)
      std::cerr << std::endl;
    //    if(debug)
    //      std::cout << std::endl;
    return sqrt(varphi/localArea);
  }

  virtual double compute_sig_ij(TriangleMesh& mesh, Vertex_point2_map &uvmap, vertex_descriptor v_i, vertex_descriptor v_j) {
    return (get(vL2Map,v_i)+get(vL2Map,v_j))/2.0;
  }

  virtual NT compute_borderLength_3D(TriangleMesh& mesh)  {
    return 1.0;
  }

};


} // namespace Surface_mesh_parameterization

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_SURFACE_MESH_PARAMETERIZATION_DISCRETE_AUTHALIC_PARAMETERIZER_3_H
