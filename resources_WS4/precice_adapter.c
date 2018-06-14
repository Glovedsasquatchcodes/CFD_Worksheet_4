#include "precice_adapter.h"
#include "boundary_val.h"

int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG){
	int dimension   = precice.getDimension();
    int getMeshID   = precice.getMeshID("FluidMesh");
    int* vertexIDs  = new int[num_coupling_cells];


    /* -------------------------CASE 1 (APPROACH 1) BEGINS------------------------------ */
	double* vertices = new double[num_coupling_cells*dimension];

    double* coords = new double[dimension];

    int coupledcellcount = 0;
	for(int j=0; j<=jmax; j++){ //left boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (0 - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
            coords = (vertices + dimension*coupledcellcount);

			precice.setMeshVertices(meshID, num_coupling_cells, coords, (vertexIDs + coupledcellcount));
			coupledcellcount++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Right boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (imax - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
            coords = (vertices + dimension*coupledcellcount);

			precice.setMeshVertices(meshID, num_coupling_cells, coords, (vertexIDs + coupledcellcount));
			coupledcellcount++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Top boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (jmax - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
            coords = (vertices + dimension*coupledcellcount);

			precice.setMeshVertices(meshID, num_coupling_cells, coords, (vertexIDs + coupledcellcount));
			coupledcellcount++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Bottom boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (0 - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;;
            coords = (vertices + dimension*coupledcellcount);

			precice.setMeshVertices(meshID, num_coupling_cells, coords, (vertexIDs + coupledcellcount));
			coupledcellcount++;
		}
		
	}
	/* -------------------------CASE 1 (APPROACH 1) ENDS------------------------------ */



    /* -------------------------CASE 1 (APPROACH 2) BEGINS------------------------------ */
	double* vertices = new double[num_coupling_cells*dimension];

    int coupledcellcount = 0;
	for(int j=0; j<=jmax; j++){ //left boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (0 - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
            coupledcellcount++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Right boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (imax - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
			coupledcellcount++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Top boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (jmax - 0.5)*dy;
			vertices[dimension*coupledcellcount + 2] = 0;
			coupledcellcount++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Bottom boundary
		if(flag[i][j]&(1<<4)){
            vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
            vertices[dimension*coupledcellcount + 1] = y_origin + (0 - 0.5)*dy;;
			vertices[dimension*coupledcellcount + 2] = 0;
			coupledcellcount++;
		}
		
	}
	precice.setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
	/* -------------------------CASE 1 (APPROACH 1) ENDS------------------------------ */


    /* -------------------CASE 2 GENERALIZED SCAN (APPROACH 1) BEGINS--------------------- */
    double* vertices = new double[num_coupling_cells*dimension];

    double* coords = new double[dimension];

    int coupledcellcount = 0;
	for(int i = 0; i<=imax; i++){
		for(int j=0; j<=jmax; j++){ 
			
			/* scanning from bottom to top, left to right */
			if(flag[i][j]&(1<<9 && (B_N(flag[i][j]) | B_S(flag[i][j])))){
				vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
				vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
				vertices[dimension*coupledcellcount + 2] = 0;
				coords = (vertices + dimension*coupledcellcount);

				precice.setMeshVertices(meshID, num_coupling_cells, coords, (vertexIDs + coupledcellcount));
				coupledcellcount++;
			}

		}
	}
    /* -------------------CASE 2 GENERALIZED SCAN (APPROACH 1) BEGINS--------------------- */


	/* -------------------CASE 2 GENERALIZED SCAN (APPROACH 2) BEGINS--------------------- */
    double* vertices = new double[num_coupling_cells*dimension];

    int coupledcellcount = 0;
	for(int i = 0; i<=imax; i++){
		for(int j=0; j<=jmax; j++){ 
			
			/* scanning from bottom to top, left to right */
			if(flag[i][j]&(1<<9 && (B_N(flag[i][j]) | B_S(flag[i][j])))){
				vertices[dimension*coupledcellcount]     = x_origin + (i - 0.5)*dx;
				vertices[dimension*coupledcellcount + 1] = y_origin + (j - 0.5)*dy;
				vertices[dimension*coupledcellcount + 2] = 0;

				coupledcellcount++;
			}

		}
	}
	precice.setMeshVertices(meshID, num_coupling_cells, vertices, vertexIDs);
				
    /* -------------------CASE 2 GENERALIZED SCAN (APPROACH 2) BEGINS--------------------- */
	free(coords);
	return vertexIDs;
	
}

void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature, int *vertexIDs,
                               int temperatureID, double **TEMP, int **FLAG)
{
	count = 0;
	for(int j=0; j<=jmax; j++){ //left boundary
		if(flag[i][j]&(1<<4)){
			temperature[count] = TEMP[1][j];
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Right boundary
		if(flag[i][j]&(1<<4)){
			temperature[count] = TEMP[imax][j];
			count++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Top boundary
		if(flag[i][j]&(1<<4)){
			temperature[count] = TEMP[i][jmax];
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Bottom boundary
		if(flag[i][j]&(1<<4)){
			temperature[count] = TEMP[i][1];
			count++;
		}
		
	}

	//precice.writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);  shift this to main
}

void set_coupling_boundary(int imax, int jmax, double dx, double dy, double *heatflux, double **TEMP, int **FLAG)
{
	//precice.readBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, heatflux); shift this to main as well
	count = 0;
	for(int j=0; j<=jmax; j++){ //left boundary
		if(flag[i][j]&(1<<4)){
			TEMP[1][j]= TEMP[0][j]+ dx*(heatflux[count]);
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Right boundary
		if(flag[i][j]&(1<<4)){
			TEMP[imax][j]= TEMP[imax+1][j]+ dx*(heatflux[count]);
			count++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Top boundary
		if(flag[i][j]&(1<<4)){
			TEMP[i][jmax]= TEMP[i][jmax+1]+ dy*(heatflux[count]);
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Bottom boundary
		if(flag[i][j]&(1<<4)){
			TEMP[i][1]= TEMP[i][0]+ dy*(heatflux[count]);
			count++;
		}
		
	}
}

	
