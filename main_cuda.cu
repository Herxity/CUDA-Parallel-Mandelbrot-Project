//Names: Ruben and Abdul
//CPU: AMD Ryzen 7 7735HS
// g++ main_threading.cpp -o ./a.out -O3 -lwebp -lpng
//(O3 Optimization enabled)
/*
mandelbrot duration: 17.9092 seconds 
*/


//GPU: NVIDIA RTX 4060 LAPTOP GPU
// nvcc main_cuda.cu -o ./a.out -O3 -lwebp -lpng
/*
mandelbrot duration: 1.61442 seconds
*/


#include <cstdint>
#include <cmath>
#include <complex>
#include <iostream>
#include <fstream>
#include <webp/encode.h>
#include <png.h>
#include <chrono>
#include <thread>
#include <vector> 
#include <cuda_runtime.h>



/*
    C = (0,0)
    f(z) = z^2 + C 
            (0,0) + (0,0)

    C = (2,0)
    f(z) = z^2 + C
            (4,0) + (2,0) = (6,0)
            (36,0) + (6,0)...


            with vecorized coding, if statements are a problem

            because all numbers in the vector are processed the same
*/
using namespace std;



//Use the __managed__ keyword to declare this array as being shared between host and device
__managed__ uint32_t* count_arr;


//Cuda does not support complex number handling
__device__ float complex_abs(float x, float y) {
    return sqrtf(x * x + y * y);
}

__device__ float mag_sqr(float x, float y) {
    return x * x + y * y;
}


__global__ void mandelbrot_worker(uint32_t* count_arr, uint32_t w, uint32_t h,
                                  const uint32_t max_count, const float xmin,
                                  const float xmax, const float ymin, const float ymax) {
    // Calculate global thread coordinates in the 2D Mandelbrot set grid
    int j = blockIdx.x * blockDim.x + threadIdx.x;  // x-coordinate
    int i = blockIdx.y * blockDim.y + threadIdx.y;  // y-coordinate

    // if (i >= h || j >= w) {
    //     printf("ERROR: OUT OF BOUNDS");
    //     return;
    // } // Bounds check

    // Map pixel coordinate (i, j) to complex plane coordinates (x, y)
    float x0 = xmin + (xmax - xmin) * j / w;
    float y0 = ymin + (ymax - ymin) * i / h;

    // Compute Mandelbrot iteration count for point (x0, y0)
    float x = 0.0f, y = 0.0f;
    uint32_t count = 0;

    while (count < max_count && complex_abs(x, y) <= 2.0f) {
        float x_temp = x * x - y * y + x0;
        y = 2.0f * x * y + y0;
        x = x_temp;
        count++;
    }

    // Write result to the output array
    count_arr[i * w + j] = count;
}

__global__ void mandelbrot_worker2(uint32_t* count_arr, uint32_t w, uint32_t h,
                                  const uint32_t max_count, const float xmin,
                                  const float dx_per_pixel, const float ymin, const float dy_per_pixel) {
    // Calculate global thread coordinates in the 2D Mandelbrot set grid
    int j = blockIdx.x * blockDim.x + threadIdx.x;  // x-coordinate
    int i = blockIdx.y * blockDim.y + threadIdx.y;  // y-coordinate

    // if (i >= h || j >= w) {
    //     printf("ERROR: OUT OF BOUNDS");
    //     return;
    // } // Bounds check

    // Map pixel coordinate (i, j) to complex plane coordinates (x, y)
    float x0 = xmin + dx_per_pixel * j ;
    float y0 = ymin + dy_per_pixel * i ;

    // Compute Mandelbrot iteration count for point (x0, y0)
    float x = 0.0f, y = 0.0f;
    uint32_t count = 0;

    while (count < max_count && mag_sqr(x, y) <= 4.0f) {
        float x_temp = x * x - y * y + x0;
        y = 2.0f * x * y + y0;
        x = x_temp;
        count++;
    }

    // Write result to the output array
    count_arr[i * w + j] = count;
}



void mandelbrot(uint32_t count_arr[], uint32_t w, uint32_t h,const uint32_t max_count, const float xmin, const float xmax, const float ymin, const float ymax){
    // sequentially write each count to array
    uint32_t num_threads = 16;
    if(num_threads%4!=0 && num_threads!=1){
        return;
    }
    // Launch kernel
    dim3 threadsPerBlock(num_threads, num_threads);  // 16x16 block of threads
    dim3 numBlocks((w + threadsPerBlock.x - 1) / threadsPerBlock.x,
                   (h + threadsPerBlock.y - 1) / threadsPerBlock.y);


    const float dx_per_pixel = ((xmax - xmin)/ w)  ;
    const float dy_per_pixel = ((ymax - ymin)/ h)  ;
    // mandelbrot_worker<<<numBlocks, threadsPerBlock>>>(count_arr, w, h, max_count, xmin, xmax, ymin, ymax);
    mandelbrot_worker2<<<numBlocks, threadsPerBlock>>>(count_arr, w, h, max_count, xmin, dx_per_pixel, ymin, dy_per_pixel);

    cudaDeviceSynchronize();//Wait for GPU code
}

void convert_mandelbrot_count_to_rgb(uint32_t pixels[], uint32_t mandelbrot_count[], uint32_t w, uint32_t h, const uint32_t colors[], uint32_t color_count) {
    for (uint32_t y = 0; y < h; y++) {
        for (uint32_t x = 0; x < w; x++) {
            uint32_t index = y * w + x;
            uint32_t count_value = mandelbrot_count[index];

            // Normalize the Mandelbrot iteration count and map it to a color
            uint32_t color_index = count_value % color_count;  // Cyclic mapping if count > color_count
            pixels[index] = colors[color_index];
        }
    }
}

void build_color_table(uint32_t colors[], uint32_t count) {
    for (uint32_t i = 0; i < count; i++) {
        // Generate a color based on the position in the palette
        uint8_t r = (i * 5) % 255;  // Adjust values to create a gradient
        uint8_t g = (i * 7) % 255;  // Feel free to tweak the multipliers
        uint8_t b = (i * 11) % 255; // to achieve different patterns
        uint8_t a = 0xFF;           // Set transparency to opaque

        // Combine color components into a single 32-bit value
        colors[i] = (a << 24) | (r << 16) | (g << 8) | b;
    }
}

bool save_webp(const char* filename, uint32_t* pixels, uint32_t w, uint32_t h, int quality) {
    // Convert the array of pixels (in RGBA format) to a WebP-encoded buffer
    uint8_t* webp_data;
    size_t webp_size = WebPEncodeRGBA((uint8_t*)pixels, w, h, w * 4, quality, &webp_data);
    
    if (webp_size == 0) {
        std::cerr << "Error encoding WebP image!" << std::endl;
        return false; // Encoding failed
    }

    // Save the WebP-encoded buffer to a file
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        std::cerr << "Error opening file for writing!" << std::endl;
        WebPFree(webp_data); // Free the WebP data in case of error
        return false;
    }

    file.write(reinterpret_cast<const char*>(webp_data), webp_size);
    file.close();
    
    // Free the WebP buffer allocated by WebPEncodeRGBA
    WebPFree(webp_data);

    return true;
}
void write_png(const char* filename, uint32_t* pixels, int w, int h) {
    FILE *fp = fopen(filename, "wb");
    if(!fp) abort();

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png) abort();

    png_infop info = png_create_info_struct(png);
    if (!info) abort();

    if (setjmp(png_jmpbuf(png))) abort();

    png_init_io(png, fp);

    // Output is 8bit depth, RGBA format.
    png_set_IHDR(
      png,
      info,
      w, h,
      8,
      PNG_COLOR_TYPE_RGBA,
      PNG_INTERLACE_NONE,
      PNG_COMPRESSION_TYPE_DEFAULT,
      PNG_FILTER_TYPE_DEFAULT
    );
    png_write_info(png, info);

    // To write image data
    png_bytep row = (png_bytep) malloc(4 * w * sizeof(png_byte));
    for(int y = 0; y < h; y++) {
        for(int x = 0; x < w; x++) {
            uint32_t pixel = pixels[y * w + x];
            png_bytep color = &(row[x * 4]);
            color[0] = (pixel >> 16) & 0xFF; // Red
            color[1] = (pixel >> 8) & 0xFF;  // Green
            color[2] = pixel & 0xFF;         // Blue
            color[3] = (pixel >> 24) & 0xFF; // Alpha
        }
        png_write_row(png, row);
    }
    png_write_end(png, NULL);

    fclose(fp);
    png_free_data(png, info, PNG_FREE_ALL, -1);
    png_destroy_write_struct(&png, &info);
    free(row);
}


int main() {
    const long w = 3840*8;  // Width for 4K * 4
    const long h = 2160*8;  // Height for 4K * 4
    uint32_t colors[64];
    uint32_t* pixels = new uint32_t[w * h];

    // Allocate managed memory
    cudaMallocManaged(&count_arr, w * h * sizeof(uint32_t));

    // Build color table
    auto start = std::chrono::high_resolution_clock::now();
    build_color_table(colors, 64);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> build_color_table_duration = end - start;
    std::cout << "build_color_table duration: " << build_color_table_duration.count() << " seconds\n";

    // Generate Mandelbrot set counts
    start = std::chrono::high_resolution_clock::now();
    mandelbrot(count_arr, w, h, 64, -2, 0.47, -1.12, 1.12);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> mandelbrot_duration = end - start;
    std::cout << "mandelbrot duration: " << mandelbrot_duration.count() << " seconds\n";
    
    // Convert counts to RGB values
    start = std::chrono::high_resolution_clock::now();
    convert_mandelbrot_count_to_rgb(pixels, count_arr, w, h, colors, 64);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> convert_to_rgb_duration = end - start;
    std::cout << "convert_mandelbrot_count_to_rgb duration: " << convert_to_rgb_duration.count() << " seconds\n";

    // Save as WebP
    start = std::chrono::high_resolution_clock::now();
    write_png("out.png", pixels, w, h);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> write_png_duration = end - start;
    std::cout << "write_png duration: " << write_png_duration.count() << " seconds\n";
    save_webp("out.webp", pixels, w, h,100);

    cudaFree(count_arr);
    delete[] pixels;
}