
namespace oomph
{


//=========================================================================
// Unstructured refineable Triangle Mesh 
//=========================================================================
template<class ELEMENT>
 class RefineableTetgenMesh : public virtual TetgenMesh<ELEMENT>,
  public virtual RefineableMeshBase
  {
   
    public:

   /// \short Build mesh, based on a TetgenMeshClosedSurface that specifies
   /// the outer boundary of the domain and any number of internal
   /// closed curves, also specified by TriangleMeshClosedSurfaces.
   /// Also specify target area for uniform element size.
   RefineableTetgenMesh(
    TetgenMeshFacetedSurface* const &outer_boundary_pt,
    Vector<TetgenMeshFacetedSurface*>& internal_closed_surface_pt,
    const double &element_volume,
    TimeStepper* time_stepper_pt=&Mesh::Default_TimeStepper,
    const bool &use_attributes=false) :
    TetgenMesh<ELEMENT>(outer_boundary_pt,
                        internal_closed_surface_pt,
                        element_volume,
                        time_stepper_pt,
                        use_attributes)
    {
     // Initialise the data associated with adaptation
      initialise_adaptation_data();

      // Setup boundary coordinates for boundaries
      /* char filename[100]; */
      /* ofstream some_file; */
      unsigned nb=nboundary();
      for (unsigned b=0;b<nb;b++)
       {
        /* sprintf(filename,"RESLT/boundary_coord_test_from_orig%i.dat",b); */
        /* some_file.open(filename); */
        this->setup_boundary_coordinates(b); //,some_file);
        /* some_file.close(); */
       }       
      //pause("done boundary coords test from orig");
    }
    
   /// \short Build mesh from specified triangulation and
   /// associated target volumes for elements in it
   RefineableTetgenMesh(const Vector<double> &target_volume,
                        tetgenio* const &tetgen_io_pt,
                        TimeStepper* time_stepper_pt=
                        &Mesh::Default_TimeStepper,
                        const bool &use_attributes=false)  
    {
     // Initialise the data associated with adaptation
     initialise_adaptation_data();
     
     // Store Timestepper used to build elements
     this->Time_stepper_pt=time_stepper_pt;
     
     // Triangulation has been created -- remember to wipe it!
     this->Tetgenio_exists =true;
     this->Tetgenio_pt = new tetgenio;

     // Add the volume constraints to the tetgenio data object
     // which may be bad because it's actually modifying things in the base
     //mesh
     //Create a local copy
     tetgenio *tetgen_input_pt = new tetgenio;;
     this->deep_copy_of_tetgenio(tetgen_io_pt,tetgen_input_pt);
     //Add volume constraints
     tetgen_input_pt->tetrahedronvolumelist = 
      new double[tetgen_input_pt->numberoftetrahedra];
     for(int e=0;e<tetgen_input_pt->numberoftetrahedra;++e)
      {
       tetgen_input_pt->tetrahedronvolumelist[e] = target_volume[e];
      }
     
     // Input string for triangle
     std::stringstream input_string_stream;
     input_string_stream<< "Vqra"; 
     
     // Convert to a *char required by the triangulate function
     char tetswitches[100];
     sprintf(tetswitches,"%s",input_string_stream.str().c_str());
     
     // Build triangulateio refined object
     tetrahedralize(tetswitches, tetgen_input_pt, this->Tetgenio_pt);       
     // Build scaffold
     this->Tmp_mesh_pt=new TetgenScaffoldMesh(*this->Tetgenio_pt);

     // Convert mesh from scaffold to actual mesh
     this->build_from_scaffold(time_stepper_pt,use_attributes);
     
     // Kill the scaffold
     delete this->Tmp_mesh_pt;
     this->Tmp_mesh_pt=0;

     //delete the input
     delete tetgen_input_pt;

     // Setup boundary coordinates for boundaries
     /* char filename[100]; */
     /* ofstream some_file; */
     unsigned nb=nboundary();
     for (unsigned b=0;b<nb;b++)
      {
       //sprintf(filename,"RESLT/boundary_coord_test_from_io%i.dat",b);
       //some_file.open(filename);
       this->setup_boundary_coordinates(b); //,some_file);
       //some_file.close();
     }       
     //pause("done boundary coords test from io");
    }   
   
   /// Empty Destructor
   virtual ~RefineableTetgenMesh() {}
   
   /// \short Problem pointer (needed for multi-domain machinery during
   /// adaptation)
   Problem*& problem_pt(){return Problem_pt;}
   
   /// Max element size allowed during adaptation
   double& max_element_size(){return Max_element_size;}
   
   /// Min element size allowed during adaptation
   double& min_element_size(){return Min_element_size;}
   
   /// Min angle before remesh gets triggered
   double& max_permitted_edge_ratio(){return Max_permitted_edge_ratio;}
   
   /// Doc the targets for mesh adaptation
   void doc_adaptivity_targets(std::ostream &outfile)
   {
    outfile << std::endl;
    outfile << "Targets for mesh adaptation: " << std::endl;
    outfile << "---------------------------- " << std::endl;
    outfile << "Target for max. error: " << Max_permitted_error << std::endl;
    outfile << "Target for min. error: " << Min_permitted_error << std::endl;
    outfile << "Target max edge ratio: " 
            << Max_permitted_edge_ratio << std::endl;
    outfile << "Min. allowed element size: " << Min_element_size << std::endl;
    outfile << "Max. allowed element size: " << Max_element_size << std::endl;
    outfile << "Don't unrefine if less than " << Max_keep_unrefined 
            << " elements need unrefinement." << std::endl;
    outfile << std::endl;
   }
   
   
   /// Refine mesh uniformly and doc process
   void refine_uniformly(DocInfo& doc_info)
   {
    throw OomphLibError("refine_uniformly() not implemented yet",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION); 
   }    
   
   
   /// \short Unrefine mesh uniformly: Return 0 for success,
   /// 1 for failure (if unrefinement has reached the coarsest permitted
   /// level)
   unsigned unrefine_uniformly()
   {
    throw OomphLibError("unrefine_uniformly() not implemented yet",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION); 
    // dummy return
    return 0;
   }
   
   /// Adapt mesh, based on elemental error provided
   void adapt(const Vector<double>& elem_error);
   
   
   //protected:
   
   /// \short Helper function that updates the input polygon's PSLG
   /// by using the end-points of elements from FaceMesh(es) that are
   /// constructed for the boundaries associated with the segments of the
   /// polygon.
   //void update_polygon_using_face_mesh(TriangleMeshPolygon* polygon_pt);
   
   /// \short Generate a new PSLG representation of the inner hole
   /// boundaries
   //virtual void surface_remesh_for_inner_hole_boundaries(
   // Vector<Vector<double> > &internal_point_coord);
   
   
   /// \short Generate a new PSLG representation of the outer boundary
   //virtual void surface_remesh_for_outer_boundary();
   
   
  /// \short Snap the boundary nodes onto any curvilinear boundaries
  //void snap_nodes_onto_boundary(RefineableTriangleMesh<ELEMENT>* &new_mesh_pt,
  //                              const unsigned &b);

   
   /// Helper function to initialise data associated with adaptation
   void initialise_adaptation_data()
   {
    // Set max and min targets for adaptation
    this->Max_element_size=1.0;
    this->Min_element_size=0.001;
    this->Max_permitted_edge_ratio=2.0;
    
    // Initialise problem pointer
    this->Problem_pt=0;
   }
   
   /// \short Build a new tetgenio object from previous TriangulateIO
   /// based on target area for each element
   //void refine_triangulateio(tetgenio& tetgen_io, 
   //                          const Vector<double> &target_volume,
   //                          tetgenio &tetgen_refine);
   

   /// \short Compute target volume based on the element's error and the
   /// error target; return max edge ratio
   double compute_volume_target(const Vector<double> &elem_error,
                                Vector<double> &target_volume)
    {
     double max_edge_ratio=0.0;
     unsigned count_unrefined=0;
     unsigned count_refined=0;
     this->Nrefinement_overruled=0;
     
     unsigned nel=this->nelement();
     for (unsigned e=0;e<nel;e++)
      {
       // Get element
       FiniteElement* el_pt=this->finite_element_pt(e);
       
       // Calculate the volume of the element
       double volume=el_pt->size();
       
       //Find the vertex coordinates
       // (vertices are enumerated first)
       double vertex[4][3];
       for(unsigned n=0;n<4;++n)
        {
         for(unsigned i=0;i<3;++i)
          {
           vertex[n][i] = el_pt->node_pt(n)->x(i);
          }
        }
       
       //Compute the radius of the circumsphere of the tetrahedron
       //Algorithm stolen from tetgen for consistency
       DenseDoubleMatrix A(3);
       for(unsigned i=0;i<3;++i)
        {
         A(0,i) = vertex[1][i] - vertex[0][i];
         A(1,i) = vertex[2][i] - vertex[0][i];
         A(2,i) = vertex[3][i] - vertex[0][i];
        }
       
       Vector<double> rhs(3);
       // Compute the right hand side vector b (3x1).
       for(unsigned i=0;i<3;++i)
        {
         rhs[i] = 0.0;
         for(unsigned k=0;k<3;++k)
          {
           rhs[i] += A(i,k)*A(i,k);
          }
         rhs[i] *= 0.5;
        }
       
       //Solve the linear system, in which the rhs is over-written with 
       //the solution
       A.solve(rhs);
       //Calculate the circum-radius
       double circum_radius = 
        sqrt(rhs[0]*rhs[0] + rhs[1]*rhs[1] + rhs[2]*rhs[2]);
       
       //Now find the shortest edge length
       Vector<double> edge(3);
       double min_length = DBL_MAX;
       for(unsigned start=0;start<4;++start)
        {
         for(unsigned end=start+1;end<4;++end)
          {
           for(unsigned i=0;i<3;++i)
            {
             edge[i] = vertex[start][i] - vertex[end][i];
            }
           double length = 
            sqrt(edge[0]*edge[0] + edge[1]*edge[1] + edge[2]*edge[2]);
           if(length < min_length) {min_length = length;}
          }
        }
       
       //Now calculate the minimum edge ratio for this element
       double local_max_edge_ratio = circum_radius/min_length;
       if(local_max_edge_ratio > max_edge_ratio)
        {
         max_edge_ratio = local_max_edge_ratio;
        }
       
       // Mimick refinement in tree-based procedure: Target volumes
       // for elements that exceed permitted error is 1/4 of their
       // current volume, corresponding to a uniform sub-division.
       if (elem_error[e] > this->max_permitted_error())
        {
         target_volume[e]=std::max(volume/4.0,Min_element_size);
         if (target_volume[e]!=Min_element_size)
          {
           count_refined++;
          }
         else
          {
           this->Nrefinement_overruled++;
          }
        }
       else if (elem_error[e] < this->min_permitted_error())
        {
         target_volume[e]=std::min(4.0*volume,Max_element_size);
         if (target_volume[e]!=Max_element_size)
          {
           count_unrefined++;
          }
        }
       else
        {
         // Leave it alone
         target_volume[e] = std::max(volume,Min_element_size); 
        }
       
      } //End of loop over elements
     
     // Tell everybody
     this->Nrefined=count_refined;
     this->Nunrefined=count_unrefined;
     
     return max_edge_ratio;
   }

   
   /// \short Problem pointer (needed for multi-domain machinery during
   /// adaptation
   Problem* Problem_pt;
   
   /// Max permitted element size
   double Max_element_size;
   
   /// Min permitted element size
   double Min_element_size;
   
   /// Max edge ratio before remesh gets triggered
   double Max_permitted_edge_ratio;
   
  }; 

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

//======================================================================
/// Adapt problem based on specified elemental error estimates
//======================================================================
template <class ELEMENT>
void RefineableTetgenMesh<ELEMENT>::adapt(const Vector<double>& elem_error)
 {    
  // Get refinement targets
  Vector<double> target_size(elem_error.size());
  double max_edge_ratio=compute_volume_target(elem_error,
                                              target_size);
  // Get maximum target volume
  unsigned n=target_size.size();
  double max_size=0.0;
  double min_size=DBL_MAX;
  for (unsigned e=0;e<n;e++)
   {
    if (target_size[e]>max_size) max_size=target_size[e];
    if (target_size[e]<min_size) min_size=target_size[e];
   }
  
  oomph_info << "Maximum target size: " << max_size << std::endl;
  oomph_info << "Minimum target size: " << min_size << std::endl;
  oomph_info << "Number of elements to be refined " 
             << this->Nrefined << std::endl;
  oomph_info << "Number of elements to be unrefined "
             << this->Nunrefined << std::endl;
  oomph_info << "Max edge ratio "<< max_edge_ratio << std::endl;

  double orig_max_size, orig_min_size;
  this->max_and_min_element_size(orig_max_size, orig_min_size);
  oomph_info << "Max/min element size in original mesh: " 
             << orig_max_size  << " "
             << orig_min_size << std::endl;    

  // Should we bother to adapt?
  if ( (Nrefined > 0) || (Nunrefined > this->max_keep_unrefined()) ||
       (max_edge_ratio > this->max_permitted_edge_ratio()) )
   {

    if (! ( (Nrefined > 0) || (Nunrefined > max_keep_unrefined()) ) )
     {
      oomph_info 
       << "Mesh regeneration triggered by edge ratio criterion\n";
     }

    //Generate a new 1D mesh representation of the inner hole boundaries
    //unsigned nhole=this->Internal_polygon_pt.size();
    //Vector<Vector<double> > internal_point_coord(nhole);
    //this->surface_remesh_for_inner_hole_boundaries(internal_point_coord);

    //Update the representation of the outer boundary
    //this->surface_remesh_for_outer_boundary();

    //If there is not a geometric object associated with the boundary
    //the reset the boundary coordinates so that the lengths are consistent
    //in the new mesh and the old mesh.
    //const  unsigned n_boundary = this->nboundary();
    //for(unsigned b=0;b<n_boundary;++b)
    // {
    //  if(this->boundary_geom_object_pt(b)==0)
    //   {
    //    this->setup_boundary_coordinates(b);
    //   }
    // }

    // Are we dealing with a solid mesh?
    SolidMesh* solid_mesh_pt=dynamic_cast<SolidMesh*>(this);

    // Build temporary uniform background mesh
    //----------------------------------------
    // with volume set by maximum required volume
    //---------------------------------------
    RefineableTetgenMesh<ELEMENT>* tmp_new_mesh_pt=0;
    /*  if (solid_mesh_pt!=0)
     {
      tmp_new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
       (closed_curve_pt,
        hole_pt,
        max_size,
        this->Time_stepper_pt,
        this->Use_attributes);
     }
     else*/
    {
     tmp_new_mesh_pt=new RefineableTetgenMesh<ELEMENT>
      (this->Outer_boundary_pt,
       this->Internal_surface_pt,
       max_size,
       this->Time_stepper_pt,
       this->Use_attributes);
    }



    // Snap to curvilinear boundaries (some code duplication as this
    // is repeated below but helper function would take so many
    // arguments that it's nearly as messy...
    
    //Pass the boundary  geometric objects to the new mesh
    //tmp_new_mesh_pt->boundary_geom_object_pt() = 
    // this->boundary_geom_object_pt();
    
    //Reset the boundary coordinates if there is
    //a geometric object associated with the boundary
    //tmp_new_mesh_pt->boundary_coordinate_limits() = 
    // this->boundary_coordinate_limits();
    //for (unsigned b=0;b<n_boundary;b++)
    // {
    //  if(tmp_new_mesh_pt->boundary_geom_object_pt(b)!=0)
    //   {
    //    tmp_new_mesh_pt->setup_boundary_coordinates(b);
    //   }
    // }
    
    //Output the mesh before any snapping takes place
    //tmp_new_mesh_pt->output("pre_mesh_nodes_snapped_0.dat");
    
    //Move the nodes on the new boundary onto the 
    //old curvilinear boundary
    //If the boundary is straight this will do precisely nothing
    //but will be somewhat inefficient
    //for(unsigned b=0;b<n_boundary;b++)
    // {
    //  this->snap_nodes_onto_boundary(tmp_new_mesh_pt,b);
    // }
    
    //Output the mesh after the snapping has taken place
    //tmp_new_mesh_pt->output("mesh_nodes_snapped_0.dat"); 
    
    // Get the tetgenio object associated with that mesh
    tetgenio *tmp_new_tetgenio_pt = tmp_new_mesh_pt->tetgenio_pt();

    
#ifdef PARANOID
    if (this->Problem_pt==0) 
     {
      throw OomphLibError("Problem pointer must be set with problem_pt()",
                          OOMPH_CURRENT_FUNCTION,
                          OOMPH_EXCEPTION_LOCATION);
     }
#endif

    RefineableTetgenMesh<ELEMENT>* new_mesh_pt=0;

    // Map storing target sizes for elements in temporary 
    // tetgenio mesh
    std::map<GeneralisedElement*,double> target_size_map;


    //////////////////////////////////////////////////////////////
    // NOTE: Repeated setup of multidomain interaction could
    // be avoided by setting up a sufficiently fine bin
    // for the original mesh and reading out the target
    // area information from there
    //////////////////////////////////////////////////////////////

    // Now start iterating to refine mesh recursively
    //-----------------------------------------------
    bool done=false;
    unsigned iter=0;
    while (!done)
     {
      
      // "Project" target volumes from current mesh onto uniform
      //------------------------------------------------------
      // background mesh
      //----------------
      
      // Temporarily switch on projection capabilities to allow
      // storage of pointer to external element.
      // Need to do this for both meshes to ensure that 
      // matching is done based on Eulerian coordinates for both
      // (in case we're dealing with solid meshes where the
      // locate_zeta would otherwise use the Lagrangian coordintes).
      unsigned nelem=this->nelement();
      for (unsigned e=0;e<nelem;e++)
       {
        dynamic_cast<ELEMENT*>(this->element_pt(e))->enable_projection();
       }
      unsigned nelem2=tmp_new_mesh_pt->nelement();
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         enable_projection();
       }

      // Set up multi domain interactions so we can figure out
      // which element in the intermediate uniform mesh is co-located
      // with given element in current mesh (which is to be refined)
      Multi_domain_functions::setup_multi_domain_interaction
       <ELEMENT>(this->Problem_pt,this,tmp_new_mesh_pt);
      
      target_size_map.clear();
      for (unsigned e=0;e<nelem;e++)
       {
        ELEMENT* el_pt=dynamic_cast<ELEMENT*>(this->element_pt(e));
        unsigned nint=el_pt->integral_pt()->nweight();
        for (unsigned ipt=0;ipt<nint;ipt++)
         {
          GeneralisedElement* ext_el_pt=el_pt->external_element_pt(0,ipt);

          // Use max. rather than min area of any element overlapping the
          // the current element, otherwise we get a rapid outward diffusion
          // of small elements
          target_size_map[ext_el_pt]=std::max(target_size_map[ext_el_pt],
                                              target_size[e]);
         }

        // Switch off projection capability          
        dynamic_cast<ELEMENT*>(this->element_pt(e))->disable_projection();
       }
      for (unsigned e=0;e<nelem2;e++)
       {
        dynamic_cast<ELEMENT*>(tmp_new_mesh_pt->element_pt(e))->
         disable_projection();
       }      

      // Now copy into target area for temporary mesh but limit to
      // the equivalent of one sub-division per iteration
      done=true;
      unsigned nel_new=tmp_new_mesh_pt->nelement();
      Vector<double> new_target_size(nel_new);
      for (unsigned e=0;e<nel_new;e++)
       {
        // No target area found for this element -- keep its size
        // by setting target area to -1 for tetrahedron
        double new_size=target_size_map[tmp_new_mesh_pt->element_pt(e)];
        if (new_size<=0.0) 
         {
          new_target_size[e]=-1.0; 
         }
        else 
         {
          // Limit target area to the equivalent of uniform
          // refinement during this stage of the iteration
          new_target_size[e]=new_size;
          if (new_target_size[e]<
              tmp_new_mesh_pt->finite_element_pt(e)->size()/4.0)
           {
            //ALH: It seems that tetgen "enlarges" the volume constraint
            //so this criterion can never be meet unless dividing by 1.2 
            //as well.
            new_target_size[e]=
             tmp_new_mesh_pt->finite_element_pt(e)->size()/4.0;
            //This is the tetgen adjustment
            new_target_size[e] /= 1.2;
         
            // We'll need to give it another go later
            done=false;
           }
         }
       }
      
      

      // Now create the new mesh from TriangulateIO structure
      //-----------------------------------------------------
      // associated with uniform background mesh and the
      //------------------------------------------------
      // associated target element sizes.
      //---------------------------------
      
      // Solid mesh?
      /*if (solid_mesh_pt!=0)
       {
        new_mesh_pt=new RefineableSolidTriangleMesh<ELEMENT>
         (new_target_area,
          tmp_new_triangulateio,
          this->Time_stepper_pt,
          this->Use_attributes);
       }      
      // No solid mesh
      else */
       { 
        new_mesh_pt=new RefineableTetgenMesh<ELEMENT>
         (new_target_size,
          tmp_new_tetgenio_pt,
          this->Time_stepper_pt,
          this->Use_attributes);
       }    
      
      
      // Snap to curvilinear boundaries (some code duplication as this
      // is repeated below but helper function would take so many
      // arguments that it's nearly as messy...
/*      
      //Pass the boundary  geometric objects to the new mesh 
      new_mesh_pt->boundary_geom_object_pt() = 
       this->boundary_geom_object_pt();
      
      
      // Reset the boundary coordinates if there is
      // a geometric object associated with the boundary
      new_mesh_pt->boundary_coordinate_limits() = 
       this->boundary_coordinate_limits();
      for (unsigned b=0;b<n_boundary;b++)
       {
        if(new_mesh_pt->boundary_geom_object_pt(b)!=0)
         {
          new_mesh_pt->setup_boundary_coordinates(b);
         }
       }
      
      //Output the mesh before any snapping takes place
      //new_mesh_pt->output("pre_mesh_nodes_snapped_1.dat"); 
      
      //Move the nodes on the new boundary onto the 
      //old curvilinear boundary
      //If the boundary is straight this will do precisely nothing
      //but will be somewhat inefficient
      for(unsigned b=0;b<n_boundary;b++)
       {
        this->snap_nodes_onto_boundary(new_mesh_pt,b);
       }
      
      //Output the mesh after the snapping has taken place
      //new_mesh_pt->output("mesh_nodes_snapped_1.dat"); 
      */
      
      // Not done: get ready for another iteration
      iter++;
      //Delete the temporary mesh
      delete tmp_new_mesh_pt;
      //Now transfer over the pointers
      if (!done)
       {
        tmp_new_mesh_pt=new_mesh_pt;
        tmp_new_tetgenio_pt=new_mesh_pt->tetgenio_pt();
       }
      
     } // end of iteration
    

    // Project current solution onto new mesh
    //---------------------------------------
    ProjectionProblem<ELEMENT>* project_problem_pt=
     new ProjectionProblem<ELEMENT>;
    project_problem_pt->mesh_pt()=new_mesh_pt;
    project_problem_pt->project(this);
    
    //this->output("pre_proj",5);
    //new_mesh_pt->output("post_proj.dat",5);
    
    //Flush the old mesh 
    unsigned nnod=nnode();
    for(unsigned j=nnod;j>0;j--)  
     { 
      delete Node_pt[j-1];  
      Node_pt[j-1] = 0; 
     } 
    unsigned nel=nelement(); 
    for(unsigned e=nel;e>0;e--)  
     { 
      delete Element_pt[e-1];  
      Element_pt[e-1] = 0; 
     } 
    
    // Now copy back to current mesh
    //------------------------------
    nnod=new_mesh_pt->nnode();
    Node_pt.resize(nnod);
    nel=new_mesh_pt->nelement();
    Element_pt.resize(nel);  
    for(unsigned j=0;j<nnod;j++)
     { 
      Node_pt[j] = new_mesh_pt->node_pt(j);
     } 
    for(unsigned e=0;e<nel;e++)
     { 
      Element_pt[e] = new_mesh_pt->element_pt(e);
     } 
    
    //Copy the boundary schemes
    unsigned nbound=new_mesh_pt->nboundary();
    Boundary_element_pt.resize(nbound);
    Face_index_at_boundary.resize(nbound);
    Boundary_node_pt.resize(nbound);
    for (unsigned b=0;b<nbound;b++)
     {
      unsigned nel=new_mesh_pt->nboundary_element(b);
      Boundary_element_pt[b].resize(nel);
      Face_index_at_boundary[b].resize(nel);
      for (unsigned e=0;e<nel;e++)
       {
        Boundary_element_pt[b][e]=new_mesh_pt->boundary_element_pt(b,e);
        Face_index_at_boundary[b][e]=new_mesh_pt->face_index_at_boundary(b,e);
       }
      unsigned nnod=new_mesh_pt->nboundary_node(b);
      Boundary_node_pt[b].resize(nnod);
      for (unsigned j=0;j<nnod;j++)
       {
        Boundary_node_pt[b][j]=new_mesh_pt->boundary_node_pt(b,j);
       }
     }

    //Also copy over the new boundary and region information
    unsigned n_region = new_mesh_pt->nregion();
    //Only bother if we have regions
    if(n_region > 1)
     {
      //Deal with the region information first
      this->Region_element_pt.resize(n_region);
      this->Region_attribute.resize(n_region);
      for(unsigned r=0;r<n_region;r++)
       {
        this->Region_attribute[r] = new_mesh_pt->region_attribute(r);
        //Find the number of elements in the region
        unsigned n_region_element = new_mesh_pt->nregion_element(r);
        this->Region_element_pt[r].resize(n_region_element);
        for(unsigned e=0;e<n_region_element;e++)
         {
          this->Region_element_pt[r][e] = new_mesh_pt->region_element_pt(r,e);
         }
       }

      //Now the boundary region information
      this->Boundary_region_element_pt.resize(nbound);
      this->Face_index_region_at_boundary.resize(nbound);
      
      //Now loop over the boundaries
      for(unsigned b=0;b<nbound;++b)
       {
        //Loop over the regions
        for(unsigned r=0;r<n_region;++r)
         {
          unsigned n_boundary_el_in_region = 
           new_mesh_pt->nboundary_element_in_region(b,r);
          
          if(n_boundary_el_in_region > 0)
           {
            //Allocate storage in the map
            this->Boundary_region_element_pt[b][r].
             resize(n_boundary_el_in_region);
            this->Face_index_region_at_boundary[b][r].
             resize(n_boundary_el_in_region);

            //Copy over the information
            for(unsigned e=0;e<n_boundary_el_in_region;++e)
             {
              this->Boundary_region_element_pt[b][r][e]
               = new_mesh_pt->boundary_element_in_region_pt(b,r,e);
              this->Face_index_region_at_boundary[b][r][e] 
               = new_mesh_pt->face_index_at_boundary_in_region(b,r,e);
             }
           }
         }
       } //End of loop over boundaries

     } //End of case when more than one region

    //Snap the newly created nodes onto any geometric objects
    //this->snap_nodes_onto_geometric_objects();

    // Copy the IDs of the vertex nodes
    //this->Oomph_vertex_nodes_id=new_mesh_pt->oomph_vertex_nodes_id();
    
    // Copy TriangulateIO representation
    //TriangleHelper::clear_triangulateio(this->Triangulateio);
    //bool quiet=true;
    //this->Triangulateio=
    // TriangleHelper::deep_copy_of_triangulateio_representation(
    //  new_mesh_pt->triangulateio_representation(),quiet);
    
    this->set_deep_copy_tetgenio_pt(new_mesh_pt->tetgenio_pt());
     //this->Tetgenio_pt = new_mesh_pt->tetgenio_pt();

    // Flush the mesh
    new_mesh_pt->flush_element_and_node_storage();
    
    // Delete the mesh and the problem
    delete new_mesh_pt;
    delete project_problem_pt;

    // Solid mesh?
    if (solid_mesh_pt!=0)
     {
      // Warning
      std::stringstream error_message;
      error_message 
       << "Lagrangian coordinates are currently not projected but are\n"
       << "are re-set during adaptation. This is not appropriate for\n"
       << "real solid mechanics problems!\n";
      OomphLibWarning(error_message.str(),
                      "RefineableTriangleMesh::adapt()",
                      OOMPH_EXCEPTION_LOCATION);
      
      // Reset Lagrangian coordinates
      dynamic_cast<SolidMesh*>(this)->set_lagrangian_nodal_coordinates();
     }
    
    double max_area;
    double min_area;
    this->max_and_min_element_size(max_area, min_area);
    oomph_info << "Max/min element size in adapted mesh: " 
               << max_area  << " "
               << min_area << std::endl;    
   }
  else
   {
    oomph_info << "Not enough benefit in adaptation.\n";
    Nrefined=0;
    Nunrefined=0;
   }
 }

}
