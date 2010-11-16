//----------------------------  extract_boundary_mesh_00.cc  ---------------------------
//    $Id: bem.cc 22693 2010-11-11 20:11:47Z kanschat $
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  extract_boundary_mesh_00.cc  ---------------------------

/*
 Code for testing the function
 GridTools::extract_boundary_mesh (...).
 We test that the order of cells and the orientation
 of the vertices is consistent between the two meshes.

 This test checks the whole thing for a 2d and 3d hypercube.
*/


#include "../tests.h"

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_out.h>
#include <grid/grid_tools.h>

using namespace std;


template <int s_dim, int spacedim>
bool test_vertices_orientation(const Triangulation<s_dim,spacedim> &boundary_mesh,
			       map< typename Triangulation<s_dim,spacedim>::cell_iterator,
				    typename Triangulation<s_dim+1,spacedim>::face_iterator >
			       &surface_to_volume_mapping,
			       const int verbosity = 1)
{
  typename Triangulation<s_dim,spacedim>::active_cell_iterator
    cell = boundary_mesh.begin_active(),
    endc = boundary_mesh.end();
  typename Triangulation<s_dim+1,spacedim>::face_iterator face;

  bool success = true;

  if (verbosity>1)
    deallog << "Vol faces" << "\t" << "Surf cell" <<
      "\t" << "Distance" <<endl;

  for (; cell!=endc; ++cell){

    face = surface_to_volume_mapping [cell];

    for (unsigned int k=0; k<GeometryInfo<s_dim>::vertices_per_cell; ++k)
      {
	Point<spacedim> diff(face->vertex(k));
	diff -= cell->vertex(k);
	if (verbosity>1){
	  deallog << face->vertex(k) << "\t\t";
	  deallog << cell->vertex(k) << "\t\t\t" << diff.square() << endl;
	}
	if (diff.square()>0){
	  success = false;
	  break;
	}
      }
    if (verbosity>1) deallog << endl;
  }
  return success;
}

template <int dim, int spacedim>
void save_mesh(const Triangulation<dim,spacedim>& tria)
{
  GridOut grid_out;
  grid_out.write_ucd (tria, deallog.get_file_stream());
}


int main ()
{

  ofstream logfile("extract_boundary_mesh_00/output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  { // Extract the whole boundary of a hyper-cube
    const int dim = 2;

    deallog << "Testing hyper_cube in dim: " << dim << "..."<< endl;
    map< Triangulation<dim-1,dim>::cell_iterator,
	 Triangulation<dim,dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_indicator(1);
    volume_mesh.refine_global (1);

    save_mesh(volume_mesh);

    Triangulation<dim-1,dim> boundary_mesh;


    GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh,
				      surface_to_volume_mapping);

    if (test_vertices_orientation(boundary_mesh, surface_to_volume_mapping))
      deallog << "Passed.";
    else
      deallog << "Failed.";
    deallog << endl;

    save_mesh(boundary_mesh);

  }

  { // Extract the whole boundary of a hyper-cube
    const int dim = 3;

    deallog << "Testing hyper_cube in dim: " << dim << "..."<< endl;
    map< Triangulation<dim-1,dim>::cell_iterator,
	 Triangulation<dim,dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_indicator(1);
    volume_mesh.refine_global (1);

    save_mesh(volume_mesh);

    Triangulation<dim-1,dim> boundary_mesh;


    GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh,
				      surface_to_volume_mapping);

    if (test_vertices_orientation(boundary_mesh, surface_to_volume_mapping))
      deallog << "Passed.";
    else
      deallog << "Failed.";
    deallog << endl;

    save_mesh(boundary_mesh);

  }


  { // Extract a piece of the boundary of a hyper-cube

    const int dim = 3;
    deallog << "Testing hyper_cube in dim: " << dim << "..."<< endl;

    map< Triangulation<dim-1,dim>::cell_iterator,
	 Triangulation<dim,dim>::face_iterator>
      surface_to_volume_mapping;

    Triangulation<dim> volume_mesh;
    GridGenerator::hyper_cube(volume_mesh);
    volume_mesh.begin_active()->face(0)->set_boundary_indicator(1);
    volume_mesh.refine_global (1);

    save_mesh(volume_mesh);

    Triangulation<dim-1,dim> boundary_mesh;
    set<unsigned char> boundary_ids;
    boundary_ids.insert(0);

    GridTools::extract_boundary_mesh (volume_mesh, boundary_mesh,
				      surface_to_volume_mapping,
				      boundary_ids);

    if (test_vertices_orientation(boundary_mesh, surface_to_volume_mapping))
      deallog << "Passed.";
    else
      deallog << "Failed.";
    deallog << endl;

    save_mesh(boundary_mesh);

  }





   return 0;
 }
