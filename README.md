# BayesfMRI
BayesfMRI R package 


### Issues to resolve:

1. Special class 
    * Make the post-processing steps and summaries separate functions to call that object (MANDY)
    * Make/edit the id_activations function to call that object (MANDY)

2. Function documentation (DAVID + MANDY)

3. Incorporate checks into functions (DAVID + MANDY)
    
4. Test on multi-session data (MANDY)

5. Rename main function and every it appears as BayesGLM


### Mesh vs. Triangles/vertices

1. Normally the user should provide the triangles and vertices
2. If the user is dealing with 2D volumetric data, they can create a mesh with inla.mesh.2d
3. If the user is dealing with 3D volumetric data, they would have to create the triangles and faces since INLA does not have a way to create a mesh in this case
4. Can we improve on the vertices provided by the user by incorporating boundaries or constraints using inla.mesh.create()?


### Masking

1. mask_vertices_faces function removes certain vertices and takes care of the faces.
2. Apparently there is a function like this in the excursions package named submesh.mesh in geometry.R (line 565).  We should look into using that instead.


### Future stuff

1. Use PARDISO library if the user has a license
2. Make the initial values for the hyperparameters and the residual precision optional
3. Translate David's create_mesh.m MATLAB code for getting vertices and faces from 3D arrays into R using deldir package
4. Consider using testthat for testing
5. Build in functionality for group analysis (RYAN?)
  
