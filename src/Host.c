#include "header.h"
#define MAX_BINARY_SIZE (0x1000000)
#define SEP Rprintf("-----------------------------------------------------------\n")
char * getKernelSource(char *filename);
const char *err_code (cl_int err_in)
{
    switch (err_in) {
        case CL_SUCCESS:
            return (char*)"CL_SUCCESS";
        case CL_DEVICE_NOT_FOUND:
            return (char*)"CL_DEVICE_NOT_FOUND";
        case CL_DEVICE_NOT_AVAILABLE:
            return (char*)"CL_DEVICE_NOT_AVAILABLE";
        case CL_COMPILER_NOT_AVAILABLE:
            return (char*)"CL_COMPILER_NOT_AVAILABLE";
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:
            return (char*)"CL_MEM_OBJECT_ALLOCATION_FAILURE";
        case CL_OUT_OF_RESOURCES:
            return (char*)"CL_OUT_OF_RESOURCES";
        case CL_OUT_OF_HOST_MEMORY:
            return (char*)"CL_OUT_OF_HOST_MEMORY";
        case CL_PROFILING_INFO_NOT_AVAILABLE:
            return (char*)"CL_PROFILING_INFO_NOT_AVAILABLE";
        case CL_MEM_COPY_OVERLAP:
            return (char*)"CL_MEM_COPY_OVERLAP";
        case CL_IMAGE_FORMAT_MISMATCH:
            return (char*)"CL_IMAGE_FORMAT_MISMATCH";
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:
            return (char*)"CL_IMAGE_FORMAT_NOT_SUPPORTED";
        case CL_BUILD_PROGRAM_FAILURE:
            return (char*)"CL_BUILD_PROGRAM_FAILURE";
        case CL_MAP_FAILURE:
            return (char*)"CL_MAP_FAILURE";
        case CL_MISALIGNED_SUB_BUFFER_OFFSET:
            return (char*)"CL_MISALIGNED_SUB_BUFFER_OFFSET";
        case CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST:
            return (char*)"CL_EXEC_STATUS_ERROR_FOR_EVENTS_IN_WAIT_LIST";
        case CL_INVALID_VALUE:
            return (char*)"CL_INVALID_VALUE";
        case CL_INVALID_DEVICE_TYPE:
            return (char*)"CL_INVALID_DEVICE_TYPE";
        case CL_INVALID_PLATFORM:
            return (char*)"CL_INVALID_PLATFORM";
        case CL_INVALID_DEVICE:
            return (char*)"CL_INVALID_DEVICE";
        case CL_INVALID_CONTEXT:
            return (char*)"CL_INVALID_CONTEXT";
        case CL_INVALID_QUEUE_PROPERTIES:
            return (char*)"CL_INVALID_QUEUE_PROPERTIES";
        case CL_INVALID_COMMAND_QUEUE:
            return (char*)"CL_INVALID_COMMAND_QUEUE";
        case CL_INVALID_HOST_PTR:
            return (char*)"CL_INVALID_HOST_PTR";
        case CL_INVALID_MEM_OBJECT:
            return (char*)"CL_INVALID_MEM_OBJECT";
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
            return (char*)"CL_INVALID_IMAGE_FORMAT_DESCRIPTOR";
        case CL_INVALID_IMAGE_SIZE:
            return (char*)"CL_INVALID_IMAGE_SIZE";
        case CL_INVALID_SAMPLER:
            return (char*)"CL_INVALID_SAMPLER";
        case CL_INVALID_BINARY:
            return (char*)"CL_INVALID_BINARY";
        case CL_INVALID_BUILD_OPTIONS:
            return (char*)"CL_INVALID_BUILD_OPTIONS";
        case CL_INVALID_PROGRAM:
            return (char*)"CL_INVALID_PROGRAM";
        case CL_INVALID_PROGRAM_EXECUTABLE:
            return (char*)"CL_INVALID_PROGRAM_EXECUTABLE";
        case CL_INVALID_KERNEL_NAME:
            return (char*)"CL_INVALID_KERNEL_NAME";
        case CL_INVALID_KERNEL_DEFINITION:
            return (char*)"CL_INVALID_KERNEL_DEFINITION";
        case CL_INVALID_KERNEL:
            return (char*)"CL_INVALID_KERNEL";
        case CL_INVALID_ARG_INDEX:
            return (char*)"CL_INVALID_ARG_INDEX";
        case CL_INVALID_ARG_VALUE:
            return (char*)"CL_INVALID_ARG_VALUE";
        case CL_INVALID_ARG_SIZE:
            return (char*)"CL_INVALID_ARG_SIZE";
        case CL_INVALID_KERNEL_ARGS:
            return (char*)"CL_INVALID_KERNEL_ARGS";
        case CL_INVALID_WORK_DIMENSION:
            return (char*)"CL_INVALID_WORK_DIMENSION";
        case CL_INVALID_WORK_GROUP_SIZE:
            return (char*)"CL_INVALID_WORK_GROUP_SIZE";
        case CL_INVALID_WORK_ITEM_SIZE:
            return (char*)"CL_INVALID_WORK_ITEM_SIZE";
        case CL_INVALID_GLOBAL_OFFSET:
            return (char*)"CL_INVALID_GLOBAL_OFFSET";
        case CL_INVALID_EVENT_WAIT_LIST:
            return (char*)"CL_INVALID_EVENT_WAIT_LIST";
        case CL_INVALID_EVENT:
            return (char*)"CL_INVALID_EVENT";
        case CL_INVALID_OPERATION:
            return (char*)"CL_INVALID_OPERATION";
        case CL_INVALID_GL_OBJECT:
            return (char*)"CL_INVALID_GL_OBJECT";
        case CL_INVALID_BUFFER_SIZE:
            return (char*)"CL_INVALID_BUFFER_SIZE";
        case CL_INVALID_MIP_LEVEL:
            return (char*)"CL_INVALID_MIP_LEVEL";
        case CL_INVALID_GLOBAL_WORK_SIZE:
            return (char*)"CL_INVALID_GLOBAL_WORK_SIZE";
        case CL_INVALID_PROPERTY:
            return (char*)"CL_INVALID_PROPERTY";
            
        default:
            return (char*)"UNKNOWN ERROR";
    }
}


void check_error(cl_int err, const char *operation, char *filename, int line)
{
    if (err != CL_SUCCESS)
    {
        Rprintf("Error during operation '%s', ", operation);
        Rprintf("in '%s' on line %d\n", filename, line);
        Rprintf("Error code was \"%s\" (%d)\n", err_code(err), err);
    }
}



unsigned getDeviceList(cl_device_id devices[MAX_DEVICES])
{
    cl_int err;
    
    // Get list of platforms
    cl_uint numPlatforms = 0;
    cl_platform_id platforms[MAX_PLATFORMS];
    err = clGetPlatformIDs(MAX_PLATFORMS, platforms, &numPlatforms);
    checkError(err, "getting platforms");
    
    // Enumerate devices
    unsigned numDevices = 0;
    for (int i = 0; i < numPlatforms; i++)
    {
        cl_uint num = 0;
        err = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL,
                             MAX_DEVICES-numDevices, devices+numDevices, &num);
        checkError(err, "getting deviceS");
        numDevices += num;
    }
    
    return numDevices;
}

void getDeviceName(cl_device_id device, char name[MAX_INFO_STRING])
{
    cl_device_info info = CL_DEVICE_NAME;
    
    // Special case for AMD
#ifdef CL_DEVICE_BOARD_NAME_AMD
    clGetDeviceInfo(device, CL_DEVICE_VENDOR, MAX_INFO_STRING, name, NULL);
    if (strstr(name, "Advanced Micro Devices"))
        info = CL_DEVICE_BOARD_NAME_AMD;
#endif
    
    clGetDeviceInfo(device, info, MAX_INFO_STRING, name, NULL);
}


int parseUInt(const char *str, cl_uint *output)
{
    char *next;
    *output = strtoul(str, &next, 10);
    return !strlen(next);
}



double sum_total(double *arr, int ngrid)
{
    double sol = 0.0;
    for (int i=0; i<ngrid; i++)
    {
        sol += arr[i];
    }
    return sol;
}


void param_OCL(int *cormod,int *NN,double *par,int *weigthed,double *nuis,int *int_par, double *dou_par)
{
    int_par[0] = cormod[0]; // Correlation Model
    int_par[1] = ncoord[0]; // number of total spatial coordinates
    int_par[2] = weigthed[0];
    int_par[3] = type[0]; //  type of distance
    int_par[4] = 8; // Size int params (power of 2)
    int_par[5] = NN[0];//
    
    dou_par[0] = par[0];
    dou_par[1] = par[1];
    dou_par[2] = par[2];
    dou_par[3] = par[3];
    dou_par[4] = nuis[0];
    dou_par[5] = nuis[1];
    dou_par[6] = maxdist[0];// the threshould of the spatial distances
    dou_par[7] = 16; // Size double params (power of 2)
    dou_par[8] = REARTH[0];// radius of the sphere
    dou_par[9] = nuis[2];
    dou_par[10] = nuis[3];
    //dou_par[11] = maxtime[0];// the threshould of the temporal distances below which the pairs are considered
}
void param_st_OCL(int *cormod,int *NN,double *par,int *weigthed,double *nuis,int *int_par, double *dou_par)
{
    int_par[0] = cormod[0]; // Correlation Model
    int_par[1] = ncoord[0]; // number of total spatial coordinates
    int_par[2] = weigthed[0];
    int_par[3] = type[0]; //  type of distance
    int_par[4] = 8; // Size int params (power of 2)
    int_par[5] = ntime[0];// number of times
    int_par[6] = NN[0];// number of times
    
    dou_par[0] = par[0];
    dou_par[1] = par[1];
    dou_par[2] = par[2];
    dou_par[3] = par[3];
    dou_par[4] = nuis[0];
    dou_par[5] = nuis[1];
    dou_par[6] = maxdist[0];// the threshould of the spatial distances
    dou_par[7] = 16; // Size double params (power of 2)
    dou_par[8] = REARTH[0];// radius of the sphere
    dou_par[9] = nuis[2];
    dou_par[10] = nuis[3];
    dou_par[11] = maxtime[0];// the threshould of the temporal distances below which the pairs are considered
    dou_par[12] = par[4];
    dou_par[13] = par[5];
    dou_par[14] = par[6];
}

void exec_kernel(double *h_x, double *h_y, double *h_mean, double *h_data, int *int_par,double *dou_par,
                 int *local_wi, int *dev, double *res, char *f_name)
{
    //#pragma OPENCL EXTEdynION cl_amd_fp64 : enable
    // Context, program, build:
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    //Rprintf("%c\n",f_name[4]);
    int length = int_par[1];
    // Vars for querying Device Info:
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list') ExecKernel!\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    
    
    
    FILE *fp;
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    fp = fopen(f_name, "r");
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    free(binary_buf);
    clBuildProgram(program, 1, &device, NULL, NULL, &err);
    kernel = clCreateKernel(program, f_name, &err);
    checkError(err, "Failed to clCreateKernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");
    //Rprintf("Maximum Kernel work-group size: %lu\n", work_group_size);
    
    
    // Creating buffers
    
    
    size_t g1;
    const int ll1 =local_wi[0];
    g1 = length + (ll1 - (length & (ll1-1))); // SPACE
    
    //Rprintf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local = ll1;
    size_t global = g1;
    
    size_t buff_size = sizeof(double) * g1;
    size_t int_par_buff = sizeof(int) * int_par[4];
    size_t dou_par_buff = sizeof(double) * dou_par[7];
    //size_t covar_buff = sizeof(double) * g1*7;
    
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordx");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, buff_size, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordx");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordy");
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, buff_size, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordy");
    
    
    cl_mem d_mean = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordy");
    err = clEnqueueWriteBuffer(commands, d_mean, CL_TRUE, 0, buff_size, (void*)h_mean, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordy");
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Writing buffer device data");
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, buff_size, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    double *sol;
    sol = (double*)calloc(g1, sizeof(double));
    cl_mem d_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device sol");
    err = clEnqueueWriteBuffer(commands, d_sol, CL_TRUE, 0, buff_size, (void*)sol, 0, NULL, NULL);
    checkError(err, "Writing buffer device sol");
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_x);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_y);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_mean);
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data);
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_sol);
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_int_par);
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &d_dou_par);
    checkError(err, "Setting kernel args length");
    
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global,&local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, d_sol, CL_TRUE, 0,buff_size, sol, 0, NULL, NULL);
    clFinish(commands);
    *res = sum_total(sol, length);
    //if(!R_FINITE(*res))  *res = LOW;
    // clean up inside kernels
    //Rprintf("GPU res: \t%f\n",res[0]);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_sol);
    clReleaseMemObject(d_mean);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
}


void print_binary_type(cl_program_binary_type program_binary_type) {
    switch (program_binary_type) {
        case CL_PROGRAM_BINARY_TYPE_NONE:
            printf("There is no binary associated with device.\n");
            break;
        case CL_PROGRAM_BINARY_TYPE_COMPILED_OBJECT:
            printf("A compiled binary is associated with device.\n");
            break;
        case CL_PROGRAM_BINARY_TYPE_LIBRARY:
            printf("A library binary is associated with device.\n");
            break;
        case CL_PROGRAM_BINARY_TYPE_EXECUTABLE:
            printf("An executable binary is associated with device.\n");
            break;
        default:
            printf("Unknown binary type.\n");
    }
}

void exec_kernel_st_dyn(double *h_x, double *h_y,double *h_t, double *h_mean, double *h_data, int *int_par,double *dou_par,
                    int *local_wi, int *dev, double *res, char *f_name,int *ns, int *NS)
{
    //Rprintf("exec_kernel_st_dyn\n");
    
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    
    FILE *fp;
    //char fileName[] = "./Kernel.clbin";
    size_t binary_size;
    char *binary_buf;
    cl_int binary_status;
    
    
    //#define MAX_BINARY_SIZE (0x100000)
    
    fp = fopen(f_name, "r");
    if (!fp) {
        Rprintf("Failed to load kernel.\n");
        
    }
    binary_buf = (char *)malloc(MAX_BINARY_SIZE);
    binary_size = fread(binary_buf, 1, MAX_BINARY_SIZE, fp);
    fclose(fp);
    
    
    program = clCreateProgramWithBinary(
                                        context, 1, &device, (const size_t *)&binary_size,
                                        (const unsigned char **)&binary_buf, &binary_status, &err
                                        );
    
    cl_program_binary_type program_binary_type;
    err = clGetProgramBuildInfo(program, device, CL_PROGRAM_BINARY_TYPE,
                                sizeof(program_binary_type), &program_binary_type, NULL);
    //print_binary_type(program_binary_type);
    
    
    clBuildProgram(program, 1, &device, NULL, NULL, &err);
    free(binary_buf);
    
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*50];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    
    kernel = clCreateKernel(program, f_name, &err);
    checkError(err, "Failed to clCreateKernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");

    int NTOT=(NS[ntime[0]-1]+ns[ntime[0]-1]);
    // Creating buffers
    size_t coords_buff = sizeof(double) * NTOT;
    size_t coordt_buff = sizeof(double) * int_par[5];//ns is a vector of lenght nt
    size_t data_buff = sizeof(double) * NTOT;
    size_t int_par_buff = sizeof(int) * int_par[4];
    size_t dou_par_buff = sizeof(double) * dou_par[7];
    
    size_t g1,g2;
    const int ll1 =local_wi[0];
    const int ll2 =local_wi[1];
   

        g1 = ncoord[0] + (ll1 - (ncoord[0] & (ll1-1))); // SPACE



  //  g1 = ncoord[0] + (ll1 - (ncoord[0] & (ll1-1))); // SPACE
    g2 = int_par[5] + (ll2 - (int_par[5] & (ll2-1))); //TIME

   
    
    //Rprintf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local[2] = {ll1,ll2};
    size_t global[2] = {g1,g2};
    
    int length = g1*g2;
    //Rprintf("LENGTH: %d\n",length);
    size_t length_buff = sizeof(double)* (length);
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    checkError(err, "Creating buffer device X");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, coords_buff, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device X");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, coords_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, coords_buff, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Creating buffer device Y");
    
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, data_buff, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    cl_mem d_t = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_t, CL_TRUE, 0, coordt_buff, (void*)h_t, 0, NULL, NULL);
    checkError(err, "Creating buffer device time");
    
    double *sol;
    sol= (double*)calloc(length, sizeof(double));
    cl_mem d_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, length_buff, NULL, &err);
    checkError(err, "Creating buffer device sol");
    err = clEnqueueWriteBuffer(commands, d_sol, CL_TRUE, 0, length_buff, (void*)sol, 0, NULL, NULL);
    checkError(err, "Writing buffer device sol");
    
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    cl_mem d_mean = clCreateBuffer(context, CL_MEM_READ_ONLY, data_buff, NULL, &err);
    checkError(err, "Creating buffer device coordy");
    err = clEnqueueWriteBuffer(commands, d_mean, CL_TRUE, 0, data_buff, (void*)h_mean, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordy");
    
    cl_mem d_ns = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    checkError(err, "Creating buffer device ns");
    err = clEnqueueWriteBuffer(commands, d_ns, CL_TRUE, 0, coordt_buff, (void*)ns, 0, NULL, NULL);
    checkError(err, "Writing buffer device ns");
    
    cl_mem d_NS = clCreateBuffer(context, CL_MEM_READ_ONLY, coordt_buff, NULL, &err);
    checkError(err, "Creating buffer device NS");
    err = clEnqueueWriteBuffer(commands, d_NS, CL_TRUE, 0, coordt_buff, (void*)NS, 0, NULL, NULL);
    checkError(err, "Writing buffer device NS");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_t); //coordt
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_x); //coordx
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_y); //coordx
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data); //data
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_mean); //data
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_sol); //mom_cond0
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &d_int_par); //mom_cond3
    err |= clSetKernelArg(kernel, 7, sizeof(cl_mem), &d_dou_par); //mom_cond3
    err |= clSetKernelArg(kernel, 8, sizeof(cl_mem), &d_ns); //ns
    err |= clSetKernelArg(kernel, 9, sizeof(cl_mem), &d_NS); //NS
    checkError(err, "Setting kernel args length");
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 2, NULL, global,local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, d_sol, CL_TRUE, 0,length_buff, sol, 0, NULL, NULL);
    clFinish(commands);
    
    *res = sum_total(sol, length);
    //Rprintf("final result: %.4f\t%.4f\t%.4f\t%.4f\n", m0,m1,m2,m3);
    
    
    // clean up inside kernels
    
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_t);
    clReleaseMemObject(d_sol);
    clReleaseMemObject(d_int_par);
    clReleaseMemObject(d_dou_par);
    clReleaseMemObject(d_ns);
    clReleaseMemObject(d_NS);
}

// ================================ Start  Create Binary Kernel
#define MAX_SOURCE_SIZE (0x100000)
#define BIN_PATH "Kernel.clbin"
#define SEP Rprintf("-----------------------------------------------------------\n")


void create_binary_kernel(int *dev, char **fname)
{
    //printf("ARCHIVO fname  %s\n", *fname);
    // Context, program, build:
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    //cl_kernel kernel;
    
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list') Compilation!\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    // Create a compute context
    /*cl_platform_id platform=NULL;
    cl_context_properties ctx_properties[] = {
        CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0
    };*/
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    

    //fclose(fp);*/
    char CL[5]; char f_nameCL[100];
    strcpy(CL, ".cl");
    strcpy(f_nameCL, *fname);
    strcat(f_nameCL, CL);
    //printf("ARCHIVO CL  %s\n", f_nameCL);
    
    
    char *kernelsource = getKernelSource(f_nameCL);
    
    program = clCreateProgramWithSource(context,1,(const char **) &kernelsource, NULL, &err);
    if(err!=CL_SUCCESS){Rprintf("Failed clCreateProgramWithSource\n");}
    
    err = clCompileProgram(program, 1, &device, "-I ./", 0, NULL, NULL, NULL, NULL);
    //err = clBuildProgram(program, 1, &device, NULL, NULL, NULL);
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048*200];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    //kernel = clCreateKernel(program, *fname, &err);
    
    
    //size_t work_group_size;
    //err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    //checkError(err, "Getting kernel work group info");
    //Rprintf("Recommended Local: %lu\n", work_group_size);
    
    
    FILE *f;
    char *binary;
    size_t binary_size;
    
    clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &binary_size, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARY_SIZES\n");}
    binary = malloc(binary_size);
    clGetProgramInfo(program, CL_PROGRAM_BINARIES, binary_size, &binary, NULL);
    //if(err!=CL_SUCCESS){Rprintf("Failed to get CL_PROGRAM_BINARIES\n");}
    f = fopen(*fname, "w");
    fwrite(binary, binary_size, 1, f);
    fclose(f);
    
    
    /* Finalization */
    clFlush(commands);
    clFinish(commands);
    clReleaseProgram(program);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    //clReleaseKernel(kernel);
    Rprintf("Successful Binary Compilation!!!\n");
}

// ================================ End  Create Binary Kernel



int DeviceInfo()
{
    
    int i, j;
    char* value;
    size_t valueSize;
    cl_uint platformCount;
    cl_platform_id* platforms;
    cl_uint deviceCount;
    cl_device_id* devices;
    cl_uint maxComputeUnits;
    cl_ulong long_entries;
    size_t p_size;
    cl_device_fp_config fp;
    // get all platforms
    clGetPlatformIDs(0, NULL, &platformCount);
    platforms = (cl_platform_id*) malloc(sizeof(cl_platform_id) * platformCount);
    clGetPlatformIDs(platformCount, platforms, NULL);
    cl_device_type dt;
    char vendor[1024];                //this strirng will hold a platforms vendor
    
    
    for (i = 0; i < platformCount; i++) {
        
        Rprintf("Platform:\t\t%u\n\n", i);
        clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, sizeof(vendor), vendor, NULL);
        Rprintf("\tPlatform Vendor:\t%s\n", vendor);
        
        // get all devices
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, 0, NULL, &deviceCount);
        devices = (cl_device_id*) malloc(sizeof(cl_device_id) * deviceCount);
        clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, deviceCount, devices, NULL);
        
        //int gpu[deviceCount];
        //int countergpu=0;
        SEP;
        // for each device print critical attributes
        for (j = 0; j < deviceCount; j++) {
            SEP;
            // print device name
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_NAME, valueSize, value, NULL);
            Rprintf("%d.\tCL_DEVICE_NAME\tDevice: \t%s\n", j, value);
            free(value);
            
            // print hardware device version
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_VERSION\tHardware version: \t%s\n", j, 1, value);
            free(value);
            
            // print software driver version
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DRIVER_VERSION\tSoftware version: \t%s\n", j, 2, value);
            free(value);
            
            // print c version supported by compiler for device
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, 0, NULL, &valueSize);
            value = (char*) malloc(valueSize);
            clGetDeviceInfo(devices[j], CL_DEVICE_OPENCL_C_VERSION, valueSize, value, NULL);
            Rprintf("%d.%d\tCL_DEVICE_OPENCL_C_VERSION\tOpenCL C version: \t%s\n", j, 3, value);
            free(value);
            
            // print parallel compute units
            clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS,
                            sizeof(maxComputeUnits), &maxComputeUnits, NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_COMPUTE_UNITS\tParallel compute units: \t%d\n", j, 4, maxComputeUnits);
            
            
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_SIZE\tGlobal Memory (MB):\t%llu\n",j, 5,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_GLOBAL_MEM_CACHE_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_GLOBAL_MEM_CACHE_SIZE\tGlobal Memory Cache (MB):\t%llu\n",j, 6,long_entries/1024/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_LOCAL_MEM_SIZE,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_LOCAL_MEM_SIZE\tLocal Memory (KB):\t%llu\n",j, 7,long_entries/1024);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_ulong),&long_entries,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_CLOCK_FREQUENCY\tMax clock (MHz) :\t%llu\n",j, 8,long_entries);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_GROUP_SIZE\tMax Work Group Size:\t%zu\n",j, 9,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_SIZES\tMax Work Item Size:\t%zu\n",j,10,p_size);
            clGetDeviceInfo(devices[j],CL_KERNEL_WORK_GROUP_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_KERNEL_WORK_GROUP_SIZE\tMax kernel Work group Size:\t%zu\n",j, 11,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_WORK_ITEM_DIMENSIONS\tMax Dev dim:\t%zu\n",j, 12,p_size);
            clGetDeviceInfo(devices[j],CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(size_t),&p_size,NULL);
            Rprintf("%d.%d\tCL_DEVICE_MAX_MEM_ALLOC_SIZE\tMax Buffer size (Mb):\t%zu\n",j, 13,p_size/1000000);
            clGetDeviceInfo(devices[j],CL_DEVICE_DOUBLE_FP_CONFIG,sizeof(cl_device_fp_config),&fp,NULL);
            Rprintf("%d.%d\tSupports double precision floating-point? %s\n",j, 14,fp != 0 ? "yes" : "no");
            clGetDeviceInfo(devices[j],CL_DEVICE_TYPE,sizeof(cl_device_type),&dt,NULL);
            Rprintf("%d.%d\tCL_DEVICE_TYPE: %s\n",j, 15,dt & CL_DEVICE_TYPE_GPU ? "GPU" : "CPU");
            //clGetDeviceInfo(devices[j],CL_DEVICE_MAX_CONSTANT_ARGS,sizeof(size_t),&p_size,NULL);
            //Rprintf("%d.%d\tCL_DEVICE_MAX_CONSTANT_ARGS\tMax kernel args:\t%zu\n",j, 15,p_size);
            SEP;
            
        }
        SEP;
        free(devices);
        
    }
    
    free(platforms);
    return 0;
    
}





char * getKernelSource(char *filename)
{
    FILE *file = fopen(filename, "r");
    if (!file)
    {
        Rprintf("Error: Could not open kernel source file\n");
        //return (EXIT_FAILURE);
    }
    fseek(file, 0, SEEK_END);
    int len = ftell(file) + 1;
    rewind(file);
    
    char *source = (char *)calloc(sizeof(char), len);
    if (!source)
    {
        Rprintf("Error: Could not allocate memory for source string\n");
        //return (EXIT_FAILURE);
    }
    fread(source, sizeof(char), len, file);
    fclose(file);
    return source;
}


void exec_kernel_source(double *h_x, double *h_y, double *h_mean, double *h_data, int *int_par,double *dou_par,
                 int *local_wi, int *dev, double *res, char *f_name)
{
    // Context, program, build:
    cl_int err;
    cl_device_id        device;     // compute device id
    cl_context       context;       // compute context
    cl_command_queue commands;      // compute command queue
    cl_program       program;       // compute program
    cl_kernel        kernel;     // compute kernel
    //Rprintf("%c\n",f_name[4]);
    int length = int_par[1];
    // Vars for querying Device Info:
    
    /*char* value;
     size_t valueSize;
     cl_uint maxComputeUnits;
     cl_ulong long_entries;
     size_t p_size;*/
    
    // Set up OpenCL context, queue, kernel, etc.
    
    cl_uint deviceIndex = dev[0];
    
    // Get list of devices
    cl_device_id devices[MAX_DEVICES];
    unsigned numDevices = getDeviceList(devices);
    
    // Check device index in range
    if (deviceIndex >= numDevices)
    {
        Rprintf("Invalid device index (try '--list')\n");
        //return EXIT_FAILURE;
    }
    
    device = devices[deviceIndex];
    // Create a compute context
    context = clCreateContext(0, 1, &device, NULL, NULL, &err);
    checkError(err, "Creating context");
    
    // Create a command queue
    commands = clCreateCommandQueue(context, device, 0, &err);
    checkError(err, "Creating command queue");
    // Create the compute program from the source buffer
    
    /*char cwd[1024];
     if (getcwd(cwd, sizeof(cwd)) != NULL)
     fRprintf(stdout, "Current working dir: %s\n", cwd);
     else
     perror("getcwd() error");*/
    
    char *kernelsource = getKernelSource("Kernel.cl");
    
    program = clCreateProgramWithSource(context, 1, (const char **) & kernelsource, NULL, &err);
    checkError(err, "Creating program");
    // Build the program
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
  
    
    
    if (err != CL_SUCCESS)
    {
        size_t len;
        char buffer[2048];
        
        Rprintf("Error: Failed to build program executable!\n%s\n", err_code(err));
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
        SEP;
        Rprintf("Build Log:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_STATUS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Status:\n%s\n", buffer);
        SEP;
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_OPTIONS, 2048*sizeof(char), buffer, &len);
        Rprintf("Build Options:\n%s\n", buffer);
        SEP;
        //return EXIT_FAILURE;
    }
    // Create the compute kernel from the program
    
    kernel = clCreateKernel(program, f_name, &err);
    checkError(err, "Creating kernel");
    
    // Find kernel work-group size
    size_t work_group_size;
    err = clGetKernelWorkGroupInfo (kernel, device, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &work_group_size, NULL);
    checkError(err, "Getting kernel work group info");
    //Rprintf("Maximum Kernel work-group size: %lu\n", work_group_size);
    
    
    // Creating buffers
    
    
    size_t g1;
    const int ll1 =local_wi[0];
    g1 = length + (ll1 - (length & (ll1-1))); // SPACE
    
    //Rprintf("GLOBAL:\t%zu\t%zu\n",g1,g2);
    size_t local = ll1;
    size_t global = g1;
    
    size_t buff_size = sizeof(double) * g1;
    size_t int_par_buff = sizeof(int) * int_par[4];
    size_t dou_par_buff = sizeof(double) * dou_par[7];
    
    
    //st_alloc = wtime();
    cl_mem d_x = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordx");
    err = clEnqueueWriteBuffer(commands, d_x, CL_TRUE, 0, buff_size, (void*)h_x, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordx");
    
    
    cl_mem d_y = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordy");
    err = clEnqueueWriteBuffer(commands, d_y, CL_TRUE, 0, buff_size, (void*)h_y, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordy");
    
    
    cl_mem d_mean = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device coordy");
    err = clEnqueueWriteBuffer(commands, d_mean, CL_TRUE, 0, buff_size, (void*)h_mean, 0, NULL, NULL);
    checkError(err, "Writing buffer device coordy");
    
    cl_mem d_data = clCreateBuffer(context, CL_MEM_READ_ONLY, buff_size, NULL, &err);
    checkError(err, "Writing buffer device data");
    err = clEnqueueWriteBuffer(commands, d_data, CL_TRUE, 0, buff_size, (void*)h_data, 0, NULL, NULL);
    checkError(err, "Creating buffer device data");
    
    double *sol;
    sol = (double*)calloc(g1, sizeof(double));
    cl_mem d_sol = clCreateBuffer(context, CL_MEM_WRITE_ONLY, buff_size, NULL, &err);
    checkError(err, "Creating buffer device sol");
    err = clEnqueueWriteBuffer(commands, d_sol, CL_TRUE, 0, buff_size, (void*)sol, 0, NULL, NULL);
    checkError(err, "Writing buffer device sol");
    
    cl_mem d_int_par = clCreateBuffer(context, CL_MEM_READ_ONLY, int_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_int_par, CL_TRUE, 0, int_par_buff, (void*)int_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device int params");
    
    
    cl_mem d_dou_par = clCreateBuffer(context, CL_MEM_READ_ONLY, dou_par_buff, NULL, &err);
    err = clEnqueueWriteBuffer(commands, d_dou_par, CL_TRUE, 0, dou_par_buff, (void*)dou_par, 0, NULL, NULL);
    checkError(err, "Creating buffer device double params");
    
    
    // Push the data out to device
    clFinish(commands);
    
    err = clSetKernelArg(kernel,  0, sizeof(cl_mem), &d_x);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &d_y);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &d_mean);
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &d_data);
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &d_sol);
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &d_int_par);
    err |= clSetKernelArg(kernel, 6, sizeof(cl_mem), &d_dou_par);
    checkError(err, "Setting kernel args length");
    
    
    // Execute the kernel
    
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global,&local, 0, NULL, NULL);
    checkError(err,"clEnqueueNDRangeKernel\n");
    clFinish(commands);
    
    
    err = clEnqueueReadBuffer(commands, d_sol, CL_TRUE, 0,buff_size, sol, 0, NULL, NULL);
    clFinish(commands);
    *res = sum_total(sol, length);
    // clean up inside kernels
    
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    
    clReleaseMemObject(d_x);
    clReleaseMemObject(d_y);
    clReleaseMemObject(d_data);
    clReleaseMemObject(d_sol);
    clReleaseMemObject(d_mean);
    //clReleaseMemObject(d_int_par);
    //clReleaseMemObject(d_dou_par);
}
