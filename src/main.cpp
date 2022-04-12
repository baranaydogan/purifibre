#include "base/config.h"
#include "image/image.h"
#include "image/sf_image.h"
#include "math/sphericalFunctions.h"
#include "tractogram/tractogramReader.h"
#include "tractogram/tractogramWriter.h"
#include "tractogram/tractogram2imageMapper.h"
#include "tractogram/gridder_4segmentLength.h"
#include "tractogram/tractogram_operators.h"

namespace CMDARGS_PURIFIBRE {
    std::string             inp_tractogram_fname;
    std::string             out_tractogram_fname;
    std::string             out_fico = "";

    float                   trimFactor  = 10;
    float                   puriFactor  = 5;

    float                   voxDim = 0;
    std::tuple<float,int>   anisotropicSmoothing (std::make_tuple(0,0));
    float                   sphericalSmoothing = 15;

    int                     numberOfThreads = 0;
    bool                    verbose         = false;
    bool                    force           = false;

}

using namespace CMDARGS_PURIFIBRE;

using namespace SF;

void purify()
{

    if ((force==false) && (existsFile(out_tractogram_fname))) {
        std::cout << "Output file, " << out_tractogram_fname << ", already exists. Use --force or -f to overwrite it." << std::endl << std::flush;
        return;
    }

    if ((force==false) && (existsFile(out_tractogram_fname))) {
        std::cout << "FICO output, " << out_fico << ", already exists. Use --force or -f to overwrite it." << std::endl << std::flush;
        return;
    }

    if (voxDim<=0) voxDim=4.0f;


    // Set default numberOfThreads
    MT::maxNumberOfThreads = (numberOfThreads==0) ? MT::maxNumberOfThreads : numberOfThreads;
    numberOfThreads        = MT::maxNumberOfThreads;

    if (!verbose) QUITE = true;

    // Initialize tractogram and make copies for multithreader
    TractogramReader* tractogram = new TractogramReader[numberOfThreads]();
    tractogram[0].initReader(inp_tractogram_fname);

    int N = tractogram[0].numberOfStreamlines;

    if (N<1) {

        if (!QUITE) std::cout << "Empty tractogram" << std::endl;
        delete[] tractogram;
        return;
    }

    for (int t = 1; t < numberOfThreads; t++)
        tractogram[t].copyFrom(tractogram[0]);
    

    // Compute sTODI
    SF::init(true,17);

    SF_Image img;
    std::vector<float> bb = getTractogramBBox(tractogram);
    bb.push_back(-0.5);
    bb.push_back(int64_t(SF::sfCoords.size())-0.5);
    img.createFromBoundingBox(4,bb,voxDim,false);

    Tractogram2ImageMapper<float> gridder(tractogram,&img);
    gridder.anisotropicSmoothing(anisotropicSmoothing);
    gridder.run(processor_4segmentLength_sf<float>, outputCompiler_4segmentLength_sf<float>);

    img.smooth(sphericalSmoothing);

    
    // Compute SECO and FICO
    std::vector<std::tuple<size_t,float>> fico;
    fico.resize(N);

    auto run = [&](MTTASK task)->void {

        int len            = tractogram[task.threadId].len[task.no];
        if (len<2) {
            fico[task.no] = std::make_tuple(task.no,0.0f);
            return;
        }

        float** streamline = tractogram[task.threadId].readStreamline(task.no);
        float T[3];

        int trim = int(len*trimFactor*0.5*0.01);
        if (trim>=int(len/2-1)) trim = int(len/2)-1;
        if (trim<=0) trim  = 0;

        float minSeco   = std::numeric_limits<float>::infinity();

        for (int l=trim; l<(len-trim-1); l++) {
            vec3sub(T,streamline[l+1],streamline[l]);
            normalize(T);

            float seco = img.getSFval(streamline[l],T); // segment-to-bundle coupling (SECO)
            if (seco < minSeco)
                minSeco = seco;
            
        }
        float seco = img.getSFval(streamline[len-trim-1],T);
        if (seco < minSeco)
            minSeco = seco;

        fico[task.no] = std::make_tuple(task.no,std::log(minSeco+1));

        for (int l=0; l<len; l++)
            delete[] streamline[l];
        delete[] streamline;

    };
    if (!QUITE) 
        MT::MTRUN(N, "Computing FICO", run);
    else
        MT::MTRUN(N, run);

    SF::clean();
    
    // Save FICO
    if (out_fico!="") {
        writeTractogram(out_fico,tractogram);

        FILE *out;
	    out = fopen(out_fico.c_str(),"ab");

        char buffer[256];

        sprintf(buffer,"CELL_DATA %lu\n",N); 	
        fwrite(buffer, sizeof(char), strlen(buffer), out);

        sprintf(buffer,"SCALARS FICO float 1\n");
        fwrite(buffer, sizeof(char), strlen(buffer), out);
        
        sprintf(buffer,"LOOKUP_TABLE default\n"); 				
        fwrite(buffer, sizeof(char), strlen(buffer), out);

        for (size_t n=0; n<N; n++) {
            float tmp = std::get<1>(fico[n]);
            swapByteOrder(tmp);
            fwrite(&tmp, sizeof(float), 1, out);
        }

        fclose(out);
    }

    for (int t = 0; t < numberOfThreads; t++) 
        tractogram[t].destroyCopy();
    delete[] tractogram;

    // Remove remN amount of smallest ones
    std::sort(fico.begin(), fico.end(), [](auto a, auto b) {return std::get<1>(a) < std::get<1>(b);} );

    std::vector<size_t> idx;
    int remN = std::floor(float(fico.size())*puriFactor*0.01);

    std::vector<std::tuple<size_t,float>> rmvdIndx;
    for (int n=0; n<remN; n++)
        rmvdIndx.push_back(fico[n]);
    std::sort(rmvdIndx.begin(), rmvdIndx.end(), [](auto a, auto b) {return std::get<0>(a) < std::get<0>(b);} );

    for (size_t n=remN; n<N; n++)
        idx.push_back(std::get<0>(fico[n]));

    // Write feature values
    writeTractogram(out_tractogram_fname,inp_tractogram_fname,idx);
    
}

int main(int argc, char *argv[])
{

    INIT();
    
    CLI::App app{"\npurifibre removes spurious streamlines from tractograms\n"};
    app.footer(FOOTER);
    app.failure_message(CLI::FailureMessage::help);
    
    app.add_option("<input tractogram>", inp_tractogram_fname, "Input tractogram (.vtk,.tck)")
        ->required()
        ->check(CLI::ExistingFile);

    app.add_option("<output tractogram>", out_tractogram_fname, "Output tractogram (.vtk)")
        ->required();
    
    app.add_option("--trim,-t", trimFactor, "Trim excludes ends of streamlines from being analyzed. E.g., when trim is 10, 90\% of the streamline is analyzed. 5\% of the streamline length from each end is excluded from the computation. Default: 10.");
    
    app.add_option("--purify, -p", puriFactor, "Percentage of streamlines to remove from the tractogram. Default: 5.");

    app.add_option("--voxDim, -v", voxDim, "Isotropic voxel dimension for sTODI computation. Default: 4.");

    app.add_option("--anisotropicSmoothing", anisotropicSmoothing, "Standard deviation of the Gaussian kernel (in mm), and computation density for anisotropic smoothing (number of streamlines). E.g. when set to 2 100, smoothing is done using 100 streamlines randomly distributed around each streamline using a Gaussion distribution with standard deviation of 2 mm. Default=0 0.");

    app.add_option("--sphericalSmoothing", sphericalSmoothing, "Amount of sTODI spherical smoothing. Default: 30.");

    app.add_option("--save_fico", out_fico, "Saves a .vtk formatted copy of input tractogram with FICO values written as a field.");

    app.add_option("--numberOfThreads,-n", numberOfThreads, "Number of threads");

    app.add_flag("--verbose", verbose, "Display info on terminal");

    app.add_flag("--force,-f", force, "Force overwriting of existing output");
    
    app.callback(purify);

    CLI11_PARSE(app, argc, argv);

    return EXIT_SUCCESS;

}