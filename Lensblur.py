import imageio
import os
import numpy as np
import scipy.ndimage

# helper functions

def each_layer_same(flat_array, output_array):
    output_array[:,:,0] = flat_array
    output_array[:,:,1] = flat_array
    output_array[:,:,2] = flat_array
    return output_array

def grayscale(image):
    averaged = np.mean(image, 2)
    gray = np.zeros_like(image)
    # gray[:,:,0] = averaged
    # gray[:,:,1] = averaged
    # gray[:,:,2] = averaged
    return each_layer_same(averaged, gray)

def normalize(image):
    gray = grayscale(image)
    min_val = np.min(gray)
    max_val = np.max(gray)
    return (gray - min_val) * 255 / max_val

def maxabs(values):
    return np.max(np.abs(values - values[4]))

def true_gradient(image):
    footprint = np.ones((3,3))
    result = np.zeros_like(image)
    result[:,:,0] = scipy.ndimage.filters.generic_filter(image[:,:,0], maxabs, footprint=footprint)
    result[:,:,1] = scipy.ndimage.filters.generic_filter(image[:,:,1], maxabs, footprint=footprint)
    result[:,:,2] = scipy.ndimage.filters.generic_filter(image[:,:,2], maxabs, footprint=footprint)
    return result

def sobel_gradient(image):
    image = grayscale(image)
    result = np.zeros_like(image)
    result = scipy.ndimage.filters.sobel(image, axis=1)
    result += scipy.ndimage.filters.sobel(image, axis=0)
    return result

def prewitt_gradient(image):
    image = grayscale(image)
    # footprint = np.ones((3,3))
    result = np.zeros_like(image)
    result= scipy.ndimage.prewitt(image, axis=1)
    result += scipy.ndimage.prewitt(image, axis=0)
    # result[:,:,1] = scipy.ndimage.filters.prewitt(image[:,:,1])
    # result[:,:,2] = scipy.ndimage.filters.prewitt(image[:,:,2])
    return result

def local_mean(image, size=5):
    final = np.zeros_like(image)
    # final[:,:,0] = scipy.ndimage.filters.convolve(image[:,:,0], weights=np.full((5,5), 1.0/25))
    # final[:,:,1] = final[:,:,0]
    # final[:,:,2] = final[:,:,0]
    return each_layer_same(scipy.ndimage.filters.convolve(image[:,:,0], weights=np.full((size,size), 1.0/(size*size))), final)

def local_median(image, size=15):
    final = np.zeros_like(image)
    # final[:,:,0] = scipy.ndimage.median_filter(image[:,:,0], (10,10))
    # final[:,:,1] = final[:,:,0]
    # final[:,:,2] = final[:,:,0]
    return each_layer_same(scipy.ndimage.median_filter(image[:,:,0], (size,size)), final)

def threshold(image, percent):
    threshold = (np.max(image) - np.min(image)) * (percent / 100.) + np.min(image)
    new_image = np.zeros_like(image)
    new_image[image < threshold] = 0
    new_image[image >= threshold] = 255
    return new_image

def laplace(image):
    final = np.zeros_like(image)
    return each_layer_same(scipy.ndimage.filters.laplace(image[:,:,0]), final)

def sharpen(image):
    final = np.zeros_like(image)
    # weights1 = np.array([[-.1, -.1, -.1],
    #                     [ -.1, 1.9, -.1],
    #                     [-.1, -.1, -.1]])
    one_layer = scipy.misc.imfilter(image[:,:,0], ftype='sharpen')
    # one_layer = scipy.misc.imfilter(one_layer, ftype='sharpen')
    # one_layer = scipy.misc.imfilter(one_layer, ftype='sharpen')
    # one_layer = scipy.misc.imfilter(one_layer, ftype='sharpen')
    # one_layer = scipy.misc.imfilter(one_layer, ftype='sharpen')
    return each_layer_same(one_layer, final)
# filepath depends on who's running the code

def NAND(image1, image2):
    final = np.zeros_like(image1)
    final[ image1 == 0 ] = 255
    final[ image2 == 0 ] = 255
    # final = np.bitwise_and(image1, image2)
    return final

filepath_will = "c:/Users/Will Molter/Pictures/Geospatial vision proj1/sample_drive/cam_3/"
filepath_jacob = "/Users/jdbruce/Downloads/WQ2017/Geospatial/sample_drive/"

filepath = filepath_jacob # change this line when you change users

def run_all(sample_drive_path):
    folders = os.listdir(sample_drive_path)
    print folders
    for folder in folders:
        if folder[:3] == "cam":
            folder_mask(sample_drive_path+folder+"/", folder)


def folder_mask(input_path, output_path):
    print "Running folder: "
    print input_path
    print "Storing output in: "
    print output_path
    # num_images = 5
    # start_image = 100

    paths = [input_path + filename for filename in os.listdir(input_path)]#[start_image: start_image + num_images]]
    prev_image = imageio.imread(paths[0])

    counter = 0
    for filename in paths[0:]:
        im = imageio.imread(filename)
        if counter == 0:
            running_total = np.zeros(imageio.imread(filename).shape)
            subtract_total = np.zeros_like(running_total)
            gradient_running_total = np.zeros_like(running_total)
            running_average = np.zeros_like(running_total)
            total_approx_dev = np.zeros_like(running_total)
        running_total += im
        running_average = (im + running_average * counter) / (counter + 1)
        total_approx_dev += np.abs(im - running_average)
        gradient_running_total += prewitt_gradient(im)
        print counter
        counter += 1
        subtract_total += np.abs(im - prev_image)
        prev_image = im

    # prep output directory
    os.mkdir(output_path)
    output_path = output_path + "/"

    # get mask 1 (from grad)

    grad_norm = normalize(gradient_running_total)
    imageio.imwrite(output_path + "grad_norm.jpg", grad_norm)
    grad_norm_threshold = threshold(grad_norm, 40)
    imageio.imwrite(output_path + "grad_norm_threshold.jpg", grad_norm_threshold)
    grad_norm_threshold_median = local_median(grad_norm_threshold)
    imageio.imwrite(output_path + "grad_norm_threshold_median.jpg", grad_norm_threshold_median)

    mask_1_final = grad_norm_threshold_median

    # get mask 2 (from dev)

    dev_normalized = normalize(total_approx_dev)
    imageio.imwrite(output_path + "deviation.jpg", dev_normalized)
    dev_normalized_median =  local_median(dev_normalized)  #imageio.imread("dev_normalized_median.jpg") #
    imageio.imwrite(output_path + "dev_normalized_median.jpg", dev_normalized_median)
    dev_normalized_median_threshold = threshold(dev_normalized_median, 10)
    imageio.imwrite(output_path + "dev_normalized_median_threshold.jpg", dev_normalized_median_threshold)
    mask_2_final = dev_normalized_median_threshold

    # get final mask by NAND of two input masks

    final_mask = NAND(mask_1_final, mask_2_final)
    imageio.imwrite(output_path + "final_mask.jpg", final_mask)
    return

run_all(filepath_jacob)
