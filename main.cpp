#include <windows.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <profileapi.h>

#include "math.h"



typedef int8_t int8;
typedef int32_t int32;
typedef int64_t int64;

typedef uint8_t uint8;
typedef uint32_t uint32;
typedef uint64_t uint64;

static bool running;

static BITMAPINFO bitmap_info;
static int bytes_per_pixel = 4;
static void *bitmap_buffer;
static int bitmap_width;
static int bitmap_height;
static bool scene_rendered = false;
static int row = 0;

#include "raytracer.cpp"


static void resize_bitmap_buffer(int width, int height) {
    if(bitmap_buffer) {
        VirtualFree(bitmap_buffer, 0, MEM_RELEASE);
    }

    bitmap_width = width;
    bitmap_height = height;
    
    bitmap_info.bmiHeader.biSize = sizeof(bitmap_info.bmiHeader);
    bitmap_info.bmiHeader.biWidth = bitmap_width;
    bitmap_info.bmiHeader.biHeight = -bitmap_height;
    bitmap_info.bmiHeader.biPlanes = 1;
    bitmap_info.bmiHeader.biBitCount = 32;
    bitmap_info.bmiHeader.biCompression = BI_RGB;

    int bitmap_buffer_size = (bitmap_width*bitmap_height)*bytes_per_pixel;
    bitmap_buffer = VirtualAlloc(0, bitmap_buffer_size, MEM_COMMIT, PAGE_READWRITE);

    row = 0;
    scene_rendered = false;
}

static void update_window(HDC DeviceContext, RECT& rect) {
    int width = rect.right - rect.left;
    int height = rect.bottom - rect.top;

    StretchDIBits(DeviceContext, 0, 0, bitmap_width, bitmap_height, 0, 0, width, height, bitmap_buffer, &bitmap_info, DIB_RGB_COLORS, SRCCOPY);
}

LRESULT CALLBACK window_callback(HWND window, UINT message, WPARAM wparam, LPARAM lparam) {
    LRESULT result = 0;

    switch(message) {
        case WM_SIZE: {
            printf("window_callback WM_SIZE\n");

            RECT rect;
            GetClientRect(window, &rect);
            int width = rect.right - rect.left;
            int height = rect.bottom - rect.top;
            resize_bitmap_buffer(width, height);
        } break;

        case WM_CLOSE: {
            running = false;
        } break;

        case WM_ACTIVATEAPP: {} break;

        case WM_DESTROY: {
            running = false;
        } break;
        
        case WM_PAINT: {
            PAINTSTRUCT Paint;
            HDC DeviceContext = BeginPaint(window, &Paint);

            RECT rect;
            GetClientRect(window, &rect);

            update_window(DeviceContext, rect);
            EndPaint(window, &Paint);
        } break;

        default: {
            result = DefWindowProc(window, message, wparam, lparam);
        } break;
    }
    
    return result;
}

int CALLBACK WinMain(HINSTANCE Instance, HINSTANCE PrevInstance, LPSTR CommandLine, int ShowCode) {
    WNDCLASS window_class = {};
    window_class.lpfnWndProc = window_callback;
    window_class.hInstance = Instance;
    window_class.lpszClassName = "Raytracer";
    bool result = AttachConsole(ATTACH_PARENT_PROCESS);

    if (!RegisterClassA(&window_class)) return 0;

    HWND window = CreateWindowExA(0, window_class.lpszClassName, "Raytracer", WS_OVERLAPPEDWINDOW|WS_VISIBLE, CW_USEDEFAULT, CW_USEDEFAULT, 640 , 480, 0, 0, Instance, 0);
    if (!window)  return 0;

    running = true;

    init_scene();

    LARGE_INTEGER frequency;
    QueryPerformanceFrequency(&frequency); 
    LARGE_INTEGER last_update;
    QueryPerformanceCounter(&last_update);
    char title[1000];

    while(running) {
        MSG message;
        while(PeekMessage(&message, 0, 0, 0, PM_REMOVE)) {
            if(message.message == WM_QUIT) {
                running = false;
            }
            
            TranslateMessage(&message);
            DispatchMessageA(&message);
        }

        if (!scene_rendered) {
            LARGE_INTEGER current_time;
            QueryPerformanceCounter(&current_time);
            int64 elapsed = (current_time.QuadPart - last_update.QuadPart) / frequency.QuadPart;
            if (elapsed >= 1) {
                sprintf(title, "%.3lg MRay/s", double(ray_count) / 1000000);
                SetWindowTextA(window, title);
                last_update = current_time;
                ray_count = 0;
            }

            scene_rendered = render_scene_row(current_camera, row++);
        }

        HDC DeviceContext = GetDC(window);
        RECT rect;
        GetClientRect(window, &rect);
        update_window(DeviceContext, rect);

        ReleaseDC(window, DeviceContext);
    }
    
    return 0;
}