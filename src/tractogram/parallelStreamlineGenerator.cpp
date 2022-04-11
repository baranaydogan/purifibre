#include "tractogram_operators.h"

void getParallelStreamlines(std::vector<std::vector<std::vector<float>>>& out, TractogramReader* _tractogram, float radius, int ringCount, int pointCountPerRing)
{

    // Initialize tractogram and make copies for multithreader
    TractogramReader* tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
        tractogram[t].setThreadId(t);
    }
    
    int N   = tractogram[0].numberOfStreamlines;
    int par = ringCount*pointCountPerRing+1;
    
    // Create output
    if (N<1)
        return;
    else
        out.resize(N*par);

    // Prepare rings
    float radStep = radius/float(ringCount);
    float angStep = TWOPI/float(pointCountPerRing);
    std::vector<float> scale_N;
    std::vector<float> scale_B;
    scale_N.push_back(0.0f);
    scale_B.push_back(0.0f);
    for (float i=1; i<(ringCount+1); i++)
        for (float j=0; j<pointCountPerRing; j++) {
            scale_N.push_back(i*radStep*std::sin(angStep*j));
            scale_B.push_back(i*radStep*std::cos(angStep*j));
        }

    // Iterate throught the whole tractogram
    auto genParallel = [&](MTTASK task)->void{

        int len = tractogram[task.threadId].len[task.no];

        for (int p=0; p<par; p++) {
            out[task.no*par+p].resize(len);
            for (int l=0; l<len; l++)
                out[task.no*par+p][l].resize(3);
        }

        if (len<2)
            return;

        float** streamline = tractogram[task.threadId].readStreamline(task.no);

        // To make cylinders around the streamlines, we will compute the parallel transport frame
        float T[3], N[3], B[3], curr_T[3];
        float proj, R_axis[3], R_angle;
        float R[4][4];

        // Get a random initial PTF and handle the initial point
        vec3sub(T,streamline[1],streamline[0]);
        normalize(T);
        MT::rndm[task.threadId]->getAUnitRandomPerpVector(N,T);
        normalize(N);
        cross(B,T,N);

        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[task.no*par+p][0][i] = streamline[0][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=1; l<(len-1); l++) {

            vec3sub(curr_T,streamline[l+1],streamline[l]);
            normalize(curr_T);

            proj       = dot(T,curr_T);
            if (proj >  1) proj =  1;
            if (proj < -1) proj = -1;
            R_angle = std::acos(proj);

            curr_T[0] += EPS6;  // for numerical stability reasons
            curr_T[1] += EPS6;
            curr_T[2] += EPS6;
            
            cross(R_axis,T,curr_T);
            normalize(R_axis);

            axisangle2Rotation(R_axis, R_angle, R);
            rotate(curr_T,T,R); // curr_T is used as a temp var
            rotate(R_axis,N,R); // R_axis is used as a temp var

            T[0] = curr_T[0];
            T[1] = curr_T[1];
            T[2] = curr_T[2];
            normalize(T);

            N[0] = R_axis[0];
            N[1] = R_axis[1];
            N[2] = R_axis[2];
            normalize(N);

            cross(B,T,N);

            for (int p=0; p<par; p++)
                for (int i=0; i<3; i++)
                    out[task.no*par+p][l][i] = streamline[l][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        }

        // Repeat the last N, B for the last segment
        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[task.no*par+p][len-1][i] = streamline[len-1][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=0; l<len; l++)
            delete[] streamline[l];
        delete[] streamline;

    };

    MT::MTRUN(N, "Generating parallel streamlines", genParallel);

}


void getParallelStreamlines(std::vector<std::vector<std::vector<float>>>& out, TractogramReader* _tractogram, float sigma, int par)
{

    // Initialize tractogram and make copies for multithreader
    TractogramReader* tractogram = new TractogramReader[MT::maxNumberOfThreads]();
    for (int t = 0; t < MT::maxNumberOfThreads; t++) {
        tractogram[t].copyFrom(*_tractogram);
        tractogram[t].setThreadId(t);
    }
    
    int N   = tractogram[0].numberOfStreamlines;
    // Create output
    if (N<1)
        return;
    else
        out.resize(par*N);

    // Prepare points to track
    std::vector<float> scale_N;
    std::vector<float> scale_B;
    scale_N.push_back(0);
    scale_B.push_back(0);
    for (int n=0; n<par; n++) {
        scale_N.push_back(MT::rndm[0]->normal_m0_s1()*sigma);
        scale_B.push_back(MT::rndm[0]->normal_m0_s1()*sigma);
    }

    // Iterate throught the whole tractogram
    auto genParallel = [&](MTTASK task)->void{

        int len = tractogram[task.threadId].len[task.no];

        for (int p=0; p<par; p++) {
            out[task.no*par+p].resize(len);
            for (int l=0; l<len; l++)
                out[task.no*par+p][l].resize(3);
        }

        if (len<2)
            return;

        float** streamline = tractogram[task.threadId].readStreamline(task.no);

        // To make cylinders around the streamlines, we will compute the parallel transport frame
        float T[3], N[3], B[3], curr_T[3];
        float proj, R_axis[3], R_angle;
        float R[4][4];

        // Get a random initial PTF and handle the initial point
        vec3sub(T,streamline[1],streamline[0]);
        normalize(T);
        MT::rndm[task.threadId]->getAUnitRandomPerpVector(N,T);
        normalize(N);
        cross(B,T,N);

        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[task.no*par+p][0][i] = streamline[0][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=1; l<(len-1); l++) {

            vec3sub(curr_T,streamline[l+1],streamline[l]);
            normalize(curr_T);

            proj       = dot(T,curr_T);
            if (proj >  1) proj =  1;
            if (proj < -1) proj = -1;
            R_angle = std::acos(proj);

            curr_T[0] += EPS6;  // for numerical stability reasons
            curr_T[1] += EPS6;
            curr_T[2] += EPS6;
            
            cross(R_axis,T,curr_T);
            normalize(R_axis);

            axisangle2Rotation(R_axis, R_angle, R);
            rotate(curr_T,T,R); // curr_T is used as a temp var
            rotate(R_axis,N,R); // R_axis is used as a temp var

            T[0] = curr_T[0];
            T[1] = curr_T[1];
            T[2] = curr_T[2];
            normalize(T);

            N[0] = R_axis[0];
            N[1] = R_axis[1];
            N[2] = R_axis[2];
            normalize(N);

            cross(B,T,N);

            for (int p=0; p<par; p++)
                for (int i=0; i<3; i++)
                    out[task.no*par+p][l][i] = streamline[l][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        }

        // Repeat the last N, B for the last segment
        for (int p=0; p<par; p++)
            for (int i=0; i<3; i++)
                out[task.no*par+p][len-1][i] = streamline[len-1][i] + N[i]*scale_N[p] + + B[i]*scale_B[p];

        for (auto l=0; l<len; l++)
            delete[] streamline[l];
        delete[] streamline;

    };

    MT::MTRUN(N, "Generating parallel streamlines", genParallel);

}