#include "raylib.h"
#include <assert.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "kiss_fft.h"
#define MINIAUDIO_IMPLEMENTATION
#include "miniaudio.h"
#define N 1024
#define RINGS 16         // Latitude (North -> South)
#define SLICES 64          // Longitude (Around) - Higher = smoother circle
kiss_fft_cfg cfg;
kiss_fft_cpx cx_in[N];       // Input (Time)
kiss_fft_cpx cx_out[N];      // Output (Frequency)
float spectrum[N];

typedef struct Point3D{
    Vector3 pos;
    Color color;
} Point3D;

const float DEFORMATION_SCALE = 4.0f;
const float WAVE_SCALE = 4.0f; 
float current_volume = 0.0f;

void data_callback(ma_device* pDevice, void* pOutput, const void* pInput, ma_uint32 frameCount)
{
    if (pInput == NULL) return;
    float* input = (float*)pInput;
    // average left and right channels to mono and store in cx_in, make sure not to exceed N
    unsigned int samples_to_process = (frameCount < N) ? frameCount : N;
    memset(cx_in, 0, sizeof(cx_in));
    for (unsigned int i = 0; i < samples_to_process; i++) {
        float sample = (input[i * 2] + input[i * 2 + 1]) / 2.0f; // assuming stereo input
        cx_in[i].r = sample;
        cx_in[i].i = 0.0f;
    }
}

void drawSphereVisualizer(Point3D points[][SLICES], float radius_base, float noise[][SLICES]){
    float time = GetTime(); // Raylib function to get seconds since start
    // --- DRAWING LOOP ---
    for (int i = 0; i < RINGS; i++)
    {
        float deformation_window = 0.0f;
        for (int j = 0; j < SLICES; j++)
        {
            // Latitude (North to South) - 0 to PI
            float phi = ((float)i / RINGS + time*0.1f + noise[i][j]) * PI; // Slight rotation over time
            // Longitude (Around) - 0 to 2PI
            float theta_base = ((float)j / SLICES + time*0.1f + noise[i][j]) * 2.0f * PI;
            float theta = theta_base + phi*0.1f; // Rotate over time + twist by latitude
            int index = (j < SLICES/2) ? j % (SLICES / 2) : SLICES / 2 - j % (SLICES / 2); // Mirror for symmetry
            float deformation = spectrum[4*index % N]; // Scale deformation
            deformation_window += deformation / 8;
            if(j > 8) deformation_window -= spectrum[(4*(j-8)) % N] / 8;
            float current_radius = radius_base + current_volume;
            float wave = sinf(time * 0.1f + theta*WAVE_SCALE + phi*WAVE_SCALE); // wave effect
            float modified_phi = phi + wave*deformation_window* DEFORMATION_SCALE;
            float modified_theta = theta + wave*deformation_window* DEFORMATION_SCALE;
            float x = current_radius * sinf(modified_phi) * cosf(modified_theta);
            float y = current_radius * cosf(modified_phi);
            float z = current_radius * sinf(modified_phi) * sinf(modified_theta);
            // Hue based on position (latitude/longitude) + spectrum influence for dynamic color shift
            float base_hue = ((float)j / SLICES) * 360.0f; // Gradient around the sphere (0-360)
            float spectrum_shift = deformation * 120.0f; // Spectrum adds color variation (shifts hue)
            float hue = fmodf(base_hue + spectrum_shift + time * 20.0f, 360.0f); // Add time for animation
            Color color = ColorFromHSV(hue, 0.8f + deformation * 0.2f, 0.9f + deformation * 0.1f);
            points[i][j].pos = (Vector3){x, y, z};
            points[i][j].color = color;
        }
    }
    for (int i = 0; i < RINGS; i++) {
        for (int j = 0; j < SLICES; j++) {
            // Draw lines to next slice
            DrawLine3D(points[i][j].pos, points[i][(j + 1) % SLICES].pos, points[i][j].color);
        }
    }
}

void initNoise(float noise[][SLICES]){
    for(int i = 0; i < RINGS; i++){
        for(int j = 0; j < SLICES; j++){
            noise[i][j] = ((float)rand()/(float)(RAND_MAX)) * 0.01f - 0.005f; // Random float between -0.005 and 0.005
        }
    }
}

int main(){
    // config audio device
    ma_device_config config = ma_device_config_init(ma_device_type_capture);
    config.capture.format   = ma_format_f32;   // We want floats
    config.capture.channels = 2;               // We want stereo input
    config.sampleRate       = 44100;           
    config.dataCallback     = data_callback;   // Use our new function
    config.pUserData        = NULL;

    ma_device device;
    if (ma_device_init(NULL, &config, &device) != MA_SUCCESS) {
        return -1; // Failed to start
    }
    
    ma_device_start(&device); // Start listening!
    // screen size
    srand(time(NULL));
    const int width = 800;
    const int height = 600;
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    InitWindow(width, height, "NCS Visualizer");
    SetTargetFPS(60);
    // Define Camera
    Camera3D camera = {0};
    camera.position = (Vector3){10.0f, 10.0f, 10.0f}; // pos
    camera.target = (Vector3){0.0f, 0.0f, 0.0f};
    camera.up = (Vector3){0.0f, 1.0f, 0.0};
    camera.fovy = 45.0f;                               // Field of view
    camera.projection = CAMERA_PERSPECTIVE;
    cfg = kiss_fft_alloc(N, 0, NULL, NULL);

    float noise[RINGS][SLICES];
    Point3D points[RINGS][SLICES];
    initNoise(noise);
    float radius_base = 5.0f;
    while(!WindowShouldClose()){ //using ESC to close 
        //Update variables here
        kiss_fft(cfg, cx_in, cx_out);
        // Calculate Magnitudes for the Visuals
        for (int i = 0; i < N; i++)
        {
            float magnitude = sqrtf(cx_out[i].r * cx_out[i].r + cx_out[i].i * cx_out[i].i);
            spectrum[i] = magnitude / N; // Update the global spectrum
        }
        float rms = 0.0f;
        for(int i=0; i<N; i++){
            rms += cx_in[i].r * cx_in[i].r;
        }
        rms = sqrtf(rms/N);
        current_volume = rms * 2.0f; // Scale it up a bit

        // UpdateCamera(&camera, CAMERA_ORBITAL);
        //Draw
        BeginDrawing();
            ClearBackground(BLACK);
            BeginMode3D(camera);
                drawSphereVisualizer(points, radius_base, noise);
                // DrawGrid(10, 1.0f);
            EndMode3D();
            DrawFPS(10, 10);
        EndDrawing();
    }
    CloseWindow(); //Close window

    return 0;
}
