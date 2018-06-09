#include "precice_adapter.h"
#include "boundary_val.h"

int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells, int meshID, int **FLAG){
	double *coords = (double *) malloc(num_coupling_cells * sizeof(int));
	count = 0;
	for(int j=0; j<=jmax; j++){ //left boundary
		if(flag[i][j]&(1<<4)){
			precice.setMeshVertices(meshID, num_coupling_cells, coords, vertexIDs);
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Right boundary
		if(flag[i][j]&(1<<4)){
			precice.setMeshVertices(meshID, num_coupling_cells, coords, vertexIDs);
			count++;
		}
		
	}
	for(int i=0; i<=imax; i++){ //Top boundary
		if(flag[i][j]&(1<<4)){
			precice.setMeshVertices(meshID, num_coupling_cells, coords, vertexIDs);
			count++;
		}
		
	}
	for(int j=0; j<=jmax; j++){ //Bottom boundary
		if(flag[i][j]&(1<<4)){
			precice.setMeshVertices(meshID, num_coupling_cells, coords, vertexIDs);
			count++;
		}
		
	}
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

	
