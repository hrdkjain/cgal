#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/mesh_segmentation.h>
#include <iostream>
#include <fstream>
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor vertex_descriptor;




int main(int argc, char** argv )
{
  Surface_mesh Mesh_3D;
  if (argc==2){
    std::ifstream input(argv[1]);
    input >> Mesh_3D;
  } else {
    std::ifstream cactus("data/cactus.off");
    cactus >> Mesh_3D;
  }
  typedef Surface_mesh::Property_map<face_descriptor,double> Facet_double_map;
  Facet_double_map face_sdf_property_map = Mesh_3D.add_property_map<face_descriptor,double>("f:sdf").first;

  CGAL::sdf_values(Mesh_3D, face_sdf_property_map);

  // create a property-map for segment-ids
  typedef Surface_mesh::Property_map<face_descriptor, std::size_t> Facet_int_map;
  Facet_int_map face_segment_property_map = Mesh_3D.add_property_map<face_descriptor,std::size_t>("f:sid").first;

  // segment the mesh using default parameters for number of levels, and smoothing lambda
  // Any other scalar values can be used instead of using SDF values computed using the CGAL function
  std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(Mesh_3D, face_sdf_property_map, face_segment_property_map);
  std::cout << "number_of_segments " <<  number_of_segments << std::endl;

  // Use this face segmentation property to mark segmentation on vertices
  typedef Surface_mesh::Property_map<vertex_descriptor, std::size_t> Vertex_id_map;
  Vertex_id_map vertex_segment_property_map = Mesh_3D.add_property_map<vertex_descriptor,std::size_t>("v:sid").first;
  typedef Surface_mesh::Property_map<vertex_descriptor, double> Vertex_sdf_map;
  Vertex_sdf_map vertex_sdf_property_map = Mesh_3D.add_property_map<vertex_descriptor,double>("v:sdf").first;
  // Assign a fix arbitary value to vertex_segment_property_map, to check for if the value is assigned to it or not
  BOOST_FOREACH(vertex_descriptor vd, vertices(Mesh_3D))
    put(vertex_segment_property_map, vd, 999);

  BOOST_FOREACH(face_descriptor fd, faces(Mesh_3D)) {
    std::size_t fid = get(face_segment_property_map, fd);
    double fsdf = get(face_sdf_property_map, fd);
    std::cout << fd << " " << fid << std::endl;
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(Mesh_3D.halfedge(fd), Mesh_3D)) {
      std::size_t vid = get(vertex_segment_property_map, vd);
      if (vid == 999) {
        //meaning that the vertex_segment_property_map is not assigned to this vertex and then we can assign the face id
        put(vertex_segment_property_map, vd, fid);
        put(vertex_sdf_property_map, vd, fsdf);
      }
      else  {
        // some value is already assigned meaning that this is the case of an vertex failin on the segmentation border
        // check for the sdf value to resolve the conflict
        if(fsdf > get(vertex_sdf_property_map, vd)) {
          // meaning that this face has more sdf then the face which previously assigned the id to it
          put(vertex_segment_property_map, vd, fid);
          put(vertex_sdf_property_map, vd, fsdf);
        }
      }
    }
  }


  std::ofstream out_FS("out.off");
  out_FS << "COFF\n";
  out_FS << Mesh_3D.number_of_vertices() << " " << Mesh_3D.number_of_faces() << " 0\n";

  // Manually save the off file with vertex color
  std::map<vertex_descriptor,std::size_t> reindex;
  int n=0;
  BOOST_FOREACH(vertex_descriptor vd, vertices(Mesh_3D))  {
    out_FS << Mesh_3D.point(vd) << " " << vertexColor(get(vertex_segment_property_map,vd)) << "\n";
    reindex[vd] = n++;
  }

  std::stringstream face_SS;
  BOOST_FOREACH(face_descriptor fd, faces(Mesh_3D))  {
    out_FS << degree(fd,Mesh_3D);
    BOOST_FOREACH(vertex_descriptor vd, vertices_around_face(Mesh_3D.halfedge(fd), Mesh_3D)) {
      out_FS << " " << reindex[vd];
    }
    out_FS << "\n";
  }

  if(out_FS.good())
    out_FS.close();
  else
    return false;
  ;
/*
  typedef CGAL::Face_filtered_graph<Surface_mesh> Filtered_graph;
  //print area of each segment and then put it in a Mesh and print it in an OFF file
  Filtered_graph segment_mesh(Mesh_3D, 0, face_segment_property_map);
  for(std::size_t id = 0; id < number_of_segments; ++id)  {
    if(id > 0)
      segment_mesh.set_selected_faces(id, face_segment_property_map);
    std::cout << "Segment "<<id<<"'s area is : "<<CGAL::Polygon_mesh_processing::area(segment_mesh)<<std::endl;
    Surface_mesh out;
    CGAL::copy_face_graph(segment_mesh, out);
    std::ostringstream oss;
    oss << "Segment_" << id<<".off";
    std::ofstream os(oss.str().data());
    os<<out;
  }
  */
}

