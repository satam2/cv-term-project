from pathlib import Path
import random

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "DW.ppm"
OUTPUT_DIR = Path(__file__).resolve().parent


def next_token(data, idx):
    n = len(data)
    while idx < n:
        c = data[idx]
        if c in b" \t\r\n":
            idx += 1
            continue
        if c == ord("#"):
            while idx < n and data[idx] not in b"\r\n":
                idx += 1
            continue
        break

    start = idx
    while idx < n and data[idx] not in b" \t\r\n#":
        idx += 1

    return data[start:idx].decode("ascii"), idx


def read_ppm(path):
    data = path.read_bytes()
    idx = 0
    magic, idx = next_token(data, idx)
    if magic != "P6":
        raise ValueError(f"Unsupported PPM format {magic}; expected P6")

    width_s, idx = next_token(data, idx)
    height_s, idx = next_token(data, idx)
    maxval_s, idx = next_token(data, idx)
    width, height, maxval = int(width_s), int(height_s), int(maxval_s)

    if maxval != 255:
        raise ValueError(f"Unsupported maxval {maxval}; expected 255")

    while idx < len(data) and data[idx] in b" \t\r\n":
        idx += 1

    pixels = data[idx:]
    expected = width * height * 3
    if len(pixels) != expected:
        raise ValueError(f"Pixel data length {len(pixels)} does not match expected {expected}")

    return width, height, bytearray(pixels)


def write_ppm(path, width, height, pixels):
    header = f"P6\n{width} {height}\n255\n".encode("ascii")
    path.write_bytes(header + bytes(pixels))


def clamp(value):
    return max(0, min(255, int(round(value))))


def gaussian_blur(width, height, pixels):
    kernel = ((1, 2, 1), (2, 4, 2), (1, 2, 1))
    output = bytearray(len(pixels))

    for y in range(height):
        for x in range(width):
            base = (y * width + x) * 3
            for channel in range(3):
                total = 0
                for ky in range(-1, 2):
                    sy = min(height - 1, max(0, y + ky))
                    for kx in range(-1, 2):
                        sx = min(width - 1, max(0, x + kx))
                        weight = kernel[ky + 1][kx + 1]
                        total += pixels[(sy * width + sx) * 3 + channel] * weight
                output[base + channel] = clamp(total / 16)

    return output


def add_gaussian_noise(pixels, sigma=25):
    output = bytearray(len(pixels))
    for i, value in enumerate(pixels):
        output[i] = clamp(value + random.gauss(0, sigma))
    return output


def main():
    random.seed(136)
    width, height, pixels = read_ppm(SOURCE)

    blurred_pixels = gaussian_blur(width, height, pixels)
    blur_path = OUTPUT_DIR / "blur.ppm"
    write_ppm(blur_path, width, height, blurred_pixels)
    print(blur_path)

    blur_noise_pixels = add_gaussian_noise(blurred_pixels)
    blur_noise_path = OUTPUT_DIR / "blur_and_noise.ppm"
    write_ppm(blur_noise_path, width, height, blur_noise_pixels)
    print(blur_noise_path)


if __name__ == "__main__":
    main()
