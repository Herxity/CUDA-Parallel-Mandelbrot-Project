# Mandelbrot Set Renderer

This project is a high-performance Mandelbrot Set renderer that utilizes both **CPU** and **GPU** implementations. The GPU implementation uses CUDA for parallel computation, significantly reducing the computation time compared to the CPU version. The rendered Mandelbrot fractals are saved as both **PNG** and **WebP** images.

In this project we compare the CPU Single Thread, Multithread, and GPU Multithreaded approaches to the problem.

---

## Features
- Efficient Mandelbrot set calculation using CUDA.
- High-resolution rendering (supports resolutions up to 4K and beyond).
- Color mapping for visually appealing fractals.
- Output formats: PNG and WebP.

---

## Build Requirements
### Hardware
- **GPU**: NVIDIA GPU with CUDA support (tested on RTX 4060 Laptop GPU).
- **CPU**: Multi-core processor (tested on AMD Ryzen 7 7735HS).

### Software
- CUDA Toolkit (nvcc compiler).
- GCC/G++ compiler.
- libpng (for PNG file output).
- libwebp (for WebP file output).

---

## Compiling the Project
### CPU Version
Compile using the following command:
```bash
g++ main_threading.cpp -o ./a.out -O3 -lwebp -lpng
```

### GPU Version
Compile using the following command:
```bash
nvcc main_cuda.cu -o ./a.out -O3 -lwebp -lpng
```

---

## Usage
Run the compiled program:
```bash
./a.out
```

### Output
- The program generates two files:
  - `out.png`: Mandelbrot fractal saved as a PNG file.
  - `out.webp`: Mandelbrot fractal saved as a WebP file.
- It also outputs execution times for different steps:
  - `build_color_table`: Time to create the color palette.
  - `mandelbrot`: Time to compute the Mandelbrot set.
  - `convert_mandelbrot_count_to_rgb`: Time to convert iteration counts to RGB values.
  - `write_png`: Time to save the PNG file.

---

## Performance
### Benchmarked on:
1. **CPU**: AMD Ryzen 7 7735HS
   - Mandelbrot calculation time: ~17.9 seconds.
2. **GPU**: NVIDIA RTX 4060 Laptop GPU
   - Mandelbrot calculation time: ~1.6 seconds.

---

## File Structure
- **`main_cuda.cu`**: GPU implementation of Mandelbrot rendering.
- **`main_threading.cpp`**: CPU implementation of Mandelbrot rendering.
- **Helper Functions**:
  - `mandelbrot_worker`: Calculates Mandelbrot iteration counts (GPU kernel).
  - `convert_mandelbrot_count_to_rgb`: Converts iteration counts to RGB color values.
  - `write_png`: Saves the image as a PNG file.
  - `save_webp`: Saves the image as a WebP file.

---

## Adjusting Parameters
### Resolution
Modify `w` and `h` in the `main()` function to change the image resolution.
Example: For 4K resolution, use:
```cpp
const long w = 3840;  // Width
const long h = 2160;  // Height
```

### Mandelbrot Parameters
Modify the bounds in the `mandelbrot` function:
```cpp
mandelbrot(count_arr, w, h, 64, -2, 0.47, -1.12, 1.12);
```
- `xmin`, `xmax`, `ymin`, `ymax`: Bounds for the complex plane.
- `max_count`: Maximum number of iterations.

### Color Palette
Modify the `build_color_table` function to experiment with different color patterns.

---

## Known Issues
- **Memory Usage**: High resolutions (e.g., 4K x 4) require significant GPU memory. Adjust resolution if memory errors occur.
- **CUDA Compatibility**: Ensure that the GPU and CUDA drivers are compatible with the program.

---

## License
This project is released under the MIT License. Feel free to use, modify, and distribute it as needed.

---

## Authors
- **Ruben**
- **Abdul**

