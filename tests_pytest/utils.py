import ctypes


def is_cuda_available():
    try:
        libcudart = ctypes.CDLL("libcudart.so")
        device_count = ctypes.c_int()
        result = libcudart.cudaGetDeviceCount(ctypes.byref(device_count))
        return result == 0 and device_count.value > 0
    except OSError:
        # libcudart.so could not be loaded
        return False
