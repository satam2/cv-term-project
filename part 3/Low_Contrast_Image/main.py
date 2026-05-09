from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SOURCE = ROOT / "Gaussian_Noise_Image" / "DW_gaussian_blur_noise.ppm"
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


def low_contrast(pixels, factor=0.55):
    output = bytearray(len(pixels))
    for i, value in enumerate(pixels):
        output[i] = clamp(128 + factor * (value - 128))
    return output


def main():
    width, height, blur_noise_pixels = read_ppm(SOURCE)

    final_pixels = low_contrast(blur_noise_pixels)
    final_path = OUTPUT_DIR / "blur_noise_and_contrast.ppm"
    write_ppm(final_path, width, height, final_pixels)
    print(final_path)


if __name__ == "__main__":
    main()
