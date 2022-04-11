#include "tractogramWriter.h"
#include <cstdint>
#include <cstdio>


bool writeTractogram(std::string fname,std::vector<std::vector<std::vector<float>>>& tractogram) {
        
    // Prepare output
    FILE *out;
	out = fopen(fname.c_str(),"wb");
	if (out==NULL) {
		std::cout << "Cannot write output. Output path doesn't exist." << std::endl;
		return false;
	}


    // Prep tractogram
    size_t streamlineCount = tractogram.size();
    size_t totalPointCount = 0;

    std::vector<int> len;
    len.resize(streamlineCount);

    for (size_t i=0; i<streamlineCount; i++) {
        len[i]           = tractogram[i].size();
        totalPointCount += len[i];
    }

    // Write header
	char buffer[256];
	sprintf(buffer, "# vtk DataFile Version 3.0\n");	    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "Generated by %s\n", SGNTR.c_str());    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "BINARY\n"); 						    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "DATASET POLYDATA\n"); 				    fwrite(buffer, sizeof(char), strlen(buffer), out);
    
    // Write points
    sprintf(buffer, "POINTS %lu float\n", totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);
	for (size_t i=0; i<streamlineCount; i++) {
        for (int p=0; p<len[i]; p++) {
            float tmp;
            tmp = tractogram[i][p][0]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = tractogram[i][p][1]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = tractogram[i][p][2]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
        }
    }

    // Write lines
    sprintf(buffer, "LINES %lu %lu\n",streamlineCount,streamlineCount+totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);

	int continue_index = 0;
	for (size_t i=0; i<streamlineCount; i++) {

		int first_count = len[i];
		swapByteOrder(first_count); fwrite(&first_count, sizeof(int), 1, out);

		int first_index	= continue_index;
		int last_index	= first_index + len[i];
		for (int i=continue_index; i<last_index; i++) {
			int tmp = i;
			swapByteOrder(tmp); fwrite(&tmp, sizeof(int), 1, out);
		}

		continue_index = last_index;
	}
 
    fclose (out);
    return true;
}




bool writeTractogram(std::string out_fname,TractogramReader* tractogram,std::vector<size_t>& idx) {

    // Prepare output
    FILE *out;
	out = fopen(out_fname.c_str(),"wb");
	if (out==NULL) {
		std::cout << "Cannot write output. Output path doesn't exist." << std::endl;
		return false;
	}

    // Prep tractogram
    size_t streamlineCount = idx.size();
    size_t totalPointCount = 0;

    std::vector<int> len;
    len.reserve(streamlineCount);

    for (size_t i=0; i<streamlineCount; i++) {
        len.push_back(tractogram->len[idx[i]]);
        totalPointCount += len[i];
    }

    // Write header
	char buffer[256];
	sprintf(buffer, "# vtk DataFile Version 3.0\n");	    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "Generated by %s\n", SGNTR.c_str());    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "BINARY\n"); 						    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "DATASET POLYDATA\n"); 				    fwrite(buffer, sizeof(char), strlen(buffer), out);
    
    // Write points
    sprintf(buffer, "POINTS %lu float\n", totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);
	for (size_t i=0; i<streamlineCount; i++) {
        float** streamline = tractogram->readStreamline(idx[i]);
        for (int p=0; p<len[i]; p++) {
            float tmp;
            tmp = streamline[p][0]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = streamline[p][1]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = streamline[p][2]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            delete[] streamline[p];
        }
        delete[] streamline;
        
    }

    // Write lines
    sprintf(buffer, "LINES %lu %lu\n",streamlineCount,streamlineCount+totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);
	int continue_index = 0;
	for (size_t i=0; i<streamlineCount; i++) {
		int first_count = len[i];
		swapByteOrder(first_count); fwrite(&first_count, sizeof(int), 1, out);

		int first_index	= continue_index;
		int last_index	= first_index + len[i];
		for (int i=continue_index; i<last_index; i++) {
			int tmp = i;
			swapByteOrder(tmp); fwrite(&tmp, sizeof(int), 1, out);
		}

		continue_index = last_index;
	}
 
    fclose (out);

    return true;

}


bool writeTractogram(std::string out_fname,TractogramReader* tractogram) {

    // Prepare output
    FILE *out;
	out = fopen(out_fname.c_str(),"wb");
	if (out==NULL) {
		std::cout << "Cannot write output. Output path doesn't exist." << std::endl;
		return false;
	}

    // Prep tractogram
    size_t streamlineCount = tractogram->numberOfStreamlines;
    size_t totalPointCount = tractogram->numberOfPoints;

    // Write header
	char buffer[256];
	sprintf(buffer, "# vtk DataFile Version 3.0\n");	    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "Generated by %s\n", SGNTR.c_str());    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "BINARY\n"); 						    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "DATASET POLYDATA\n"); 				    fwrite(buffer, sizeof(char), strlen(buffer), out);
    
    // Write points
    sprintf(buffer, "POINTS %lu float\n", totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);
	for (size_t i=0; i<streamlineCount; i++) {
        float** streamline = tractogram->readStreamline(i);
        for (uint32_t p=0; p<tractogram->len[i]; p++) {
            float tmp;
            tmp = streamline[p][0]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = streamline[p][1]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = streamline[p][2]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            delete[] streamline[p];
        }
        delete[] streamline;
        
    }

    // Write lines
    sprintf(buffer, "LINES %lu %lu\n",streamlineCount,streamlineCount+totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);

	int continue_index = 0;
	for (size_t i=0; i<streamlineCount; i++) {

		int first_count = tractogram->len[i];
		swapByteOrder(first_count); fwrite(&first_count, sizeof(int), 1, out);

		int first_index	= continue_index;
		int last_index	= first_index + tractogram->len[i];
		for (int i=continue_index; i<last_index; i++) {
			int tmp = i;
			swapByteOrder(tmp); fwrite(&tmp, sizeof(int), 1, out);
		}

		continue_index = last_index;
	}
 
    fclose (out);
    return true;

}



bool writeTractogram(std::string out_fname,std::string inp_fname,std::vector<size_t>& idx) {
    TractogramReader tractogram(inp_fname);
    return writeTractogram(out_fname,&tractogram,idx);    
}


bool writeTractogram(std::string out_kept_fname,std::string out_rmvd_fname,std::string inp_fname,std::vector<size_t>& idx) {

    // Prep input
    
    TractogramReader tractogram(inp_fname);
    writeTractogram(out_kept_fname,&tractogram,idx);
    

    // Rmvd tractogram
    std::vector<size_t> rmvd_idx;
    for (size_t i=0; i<tractogram.numberOfStreamlines; i++)
        rmvd_idx.push_back(i);

    removeIdx(rmvd_idx,idx);

    writeTractogram(out_rmvd_fname,&tractogram,rmvd_idx);

    return true;
}
















bool writeTractogram(std::string fname,std::vector<std::vector<std::vector<float>>>& tractogram,std::vector<TractogramField>& fields) {
        
    // Prepare output
    FILE *out;
	out = fopen(fname.c_str(),"wb");
	if (out==NULL) {
		std::cout << "Cannot write output. Output path doesn't exist." << std::endl;
		return false;
	}


    // Prep tractogram
    size_t streamlineCount = tractogram.size();
    size_t totalPointCount = 0;

    std::vector<int> len;
    len.resize(streamlineCount);

    for (size_t i=0; i<streamlineCount; i++) {
        len[i]           = tractogram[i].size();
        totalPointCount += len[i];
    }

    // Write header
	char buffer[256];
	sprintf(buffer, "# vtk DataFile Version 3.0\n");	    fwrite(buffer, sizeof(char), strlen(buffer), out);
    sprintf(buffer, "Generated by %s\n", SGNTR.c_str());    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "BINARY\n"); 						    fwrite(buffer, sizeof(char), strlen(buffer), out);
	sprintf(buffer, "DATASET POLYDATA\n"); 				    fwrite(buffer, sizeof(char), strlen(buffer), out);
    
    // Write points
    sprintf(buffer, "POINTS %lu float\n", totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);
	for (size_t i=0; i<streamlineCount; i++) {
        for (int p=0; p<len[i]; p++) {
            float tmp;
            tmp = tractogram[i][p][0]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = tractogram[i][p][1]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
            tmp = tractogram[i][p][2]; 	swapByteOrder(tmp); fwrite(&tmp, sizeof(float), 1, out);
        }
    }

    // Write lines
    sprintf(buffer, "LINES %lu %lu\n",streamlineCount,streamlineCount+totalPointCount);
	fwrite(buffer, sizeof(char), strlen(buffer), out);

	int continue_index = 0;
	for (size_t i=0; i<streamlineCount; i++) {

		int first_count = len[i];
		swapByteOrder(first_count); fwrite(&first_count, sizeof(int), 1, out);

		int first_index	= continue_index;
		int last_index	= first_index + len[i];
		for (int i=continue_index; i<last_index; i++) {
			int tmp = i;
			swapByteOrder(tmp); fwrite(&tmp, sizeof(int), 1, out);
		}

		continue_index = last_index;
	}



    // Write fields
    std::vector<int> cellDataInd;
    std::vector<int> pointDataInd;
    
    for (size_t i=0; i<fields.size(); i++) {
        if (fields[i].owner == STREAMLINE) cellDataInd.push_back(i);
        if (fields[i].owner == POINT)      pointDataInd.push_back(i);
    }
    
    // Write STREAMLINE fields
    if (cellDataInd.size()>0) {
        sprintf(buffer,"CELL_DATA %lu\n",streamlineCount); 	fwrite(buffer, sizeof(char), strlen(buffer), out);
        
        for (size_t i=0; i<cellDataInd.size(); i++) {
            
            sprintf(buffer,"SCALARS %s ",fields[cellDataInd[i]].name.c_str());
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            if (fields[cellDataInd[i]].datatype == FLOAT32) sprintf(buffer,"float ");
            if (fields[cellDataInd[i]].datatype == INT32)   sprintf(buffer,"int ");
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            sprintf(buffer,"%d\n",fields[cellDataInd[i]].dimension); 	
            fwrite(buffer, sizeof(char), strlen(buffer), out);
            
            sprintf(buffer,"LOOKUP_TABLE default\n"); 				
            fwrite(buffer, sizeof(char), strlen(buffer), out);
            
            for (size_t s=0; s<streamlineCount; s++) {
                for (int d=0; d<fields[cellDataInd[i]].dimension; d++) {
                    
                    if (fields[cellDataInd[i]].datatype==INT32) {
                        int** data = (int**)(fields[cellDataInd[i]].data[s]);
                        int tmp = data[0][d];
                        swapByteOrder(tmp);
                        fwrite(&tmp, sizeof(int), 1, out);
                    }
                    if (fields[cellDataInd[i]].datatype==FLOAT32) {
                        float** data = (float**)(fields[cellDataInd[i]].data[s]);
                        float tmp = data[0][d];
                        swapByteOrder(tmp);
                        fwrite(&tmp, sizeof(float), 1, out);
                    }
                    
                }
            }
            
        }
        
    }
    
    // Write POINT fields
    if (pointDataInd.size()>0) {
        sprintf(buffer,"POINT_DATA %lu\n",totalPointCount); 	fwrite(buffer, sizeof(char), strlen(buffer), out);
        
        for (size_t i=0; i<pointDataInd.size(); i++) {
            
            sprintf(buffer,"SCALARS %s ",fields[cellDataInd[i]].name.c_str());
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            if (fields[pointDataInd[i]].datatype == FLOAT32) sprintf(buffer,"float ");
            if (fields[pointDataInd[i]].datatype == INT32)   sprintf(buffer,"int ");
            fwrite(buffer, sizeof(char), strlen(buffer), out);

            sprintf(buffer,"%d\n",fields[pointDataInd[i]].dimension); 	
            fwrite(buffer, sizeof(char), strlen(buffer), out);
            
            sprintf(buffer,"LOOKUP_TABLE default\n");
            fwrite(buffer, sizeof(char), strlen(buffer), out);
            
            for (size_t s=0; s<streamlineCount; s++) {
                for (int l=0; l<len[s]; l++) {
                    for (int d=0; d<fields[pointDataInd[i]].dimension; d++) {

                        if (fields[pointDataInd[i]].datatype==INT32) {
                            int** data = (int**)(fields[pointDataInd[i]].data[s]);
                            int tmp = data[l][d];
                            swapByteOrder(tmp);
                            fwrite(&tmp, sizeof(int), 1, out);
                        }

                        if (fields[pointDataInd[i]].datatype==FLOAT32) {
                            float** data = (float**)(fields[pointDataInd[i]].data[s]);
                            int tmp = data[l][d];
                            swapByteOrder(tmp);
                            fwrite(&tmp, sizeof(float), 1, out);

                        }

                    }
                }
            }
            
        }
        
    }




 
    fclose (out);
    return true;
}