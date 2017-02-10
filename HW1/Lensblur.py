import imageio
import os
import numpy as np
import scipy.ndimage
import sys

# Automatic Lens Smear Detection
# Geospatial Vision and Visualization
# WQ 2017 HW1
# Will Molter and Jacob Bruce
#
# To run from the command line, provide the path to the
# sample_drive directory (the expanded .rar) which holds
# the folders cam_0, 1, 2, 3, and 5.
# Expected running time is around 2 hours per folder
# Example: python Lensblur.py /Users/jdbruce/Downloads/WQ2017/Geospatial/sample_drive/

# helper functions

def each_layer_same(flat_array, output_array):
    # takes 2D array of values and copies it into each layer of a 3D array
    # useful for performing grayscale operations on one layer and then replacing it
    output_array[:,:,0] = flat_array
    output_array[:,:,1] = flat_array
    output_array[:,:,2] = flat_array
    return output_array

def grayscale(image):
    # gives the grayscale of the image by averaging across RGB
    averaged = np.mean(image, 2)
    gray = np.zeros_like(image)
    return each_layer_same(averaged, gray)

def normalize(image):
    # normalizes a 3D array to be written as an image
    # stretches it to put its minimum value at 0 and its maximum at xFF
    gray = grayscale(image)
    min_val = np.min(gray)
    max_val = np.max(gray)
    return (gray - min_val) * 255 / max_val

def maxabs(values):
    # gives the maximum difference between the center value and any of the others
    # useful for calculating the gradient
    return np.max(np.abs(values - values[4]))

def true_gradient(image):
    # calculates the gradient of an image by calculating the maximum change at each pixel
    footprint = np.ones((3,3))
    result = np.zeros_like(image)
    result[:,:,0] = scipy.ndimage.filters.generic_filter(image[:,:,0], maxabs, footprint=footprint)
    result[:,:,1] = scipy.ndimage.filters.generic_filter(image[:,:,1], maxabs, footprint=footprint)
    result[:,:,2] = scipy.ndimage.filters.generic_filter(image[:,:,2], maxabs, footprint=footprint)
    return result

def sobel_gradient(image):
    # wrapper function for scipy's sobel gradient to run on images
    image = grayscale(image)
    result = np.zeros_like(image)
    result = scipy.ndimage.filters.sobel(image, axis=1)
    result += scipy.ndimage.filters.sobel(image, axis=0)
    return result

def prewitt_gradient(image):
    # wrapper function for scipy's prewitt gradient to run on images
    image = grayscale(image)
    result = np.zeros_like(image)
    result= scipy.ndimage.prewitt(image, axis=1)
    result += scipy.ndimage.prewitt(image, axis=0)
    return result

def local_mean(image, size=5):
    # caclulates the local mean for a pixel within a neighborhood of a given size
    final = np.zeros_like(image)
    return each_layer_same(scipy.ndimage.filters.convolve(image[:,:,0], weights=np.full((size,size), 1.0/(size*size))), final)

def local_median(image, size=15):
    # caclulates the local median for a pixel within a neighborhood of a given size
    final = np.zeros_like(image)
    return each_layer_same(scipy.ndimage.median_filter(image[:,:,0], (size,size)), final)

def threshold(image, percent):
    # turns the image into a black and white image by setting a percent of the image range
    # above which pixels are only white and below which pixels are only black
    # should only be used on grayscale images
    threshold = (np.max(image) - np.min(image)) * (percent / 100.) + np.min(image)
    new_image = np.zeros_like(image)
    new_image[image < threshold] = 0
    new_image[image >= threshold] = 255
    return new_image

def laplace(image):
    # wrapper function for scipy's laplace operator to run on images
    final = np.zeros_like(image)
    return each_layer_same(scipy.ndimage.filters.laplace(image[:,:,0]), final)

def sharpen(image):
    # wrapper function for scipy's sharpen operator to run on images
    final = np.zeros_like(image)
    one_layer = scipy.misc.imfilter(image[:,:,0], ftype='sharpen')
    return each_layer_same(one_layer, final)

def NAND(image1, image2):
    # takes two inverted masks and returns their NAND
    # assuming all inputs are already thresholded to 0 and xFF
    final = np.zeros_like(image1)
    final[ image1 == 0 ] = 255
    final[ image2 == 0 ] = 255
    return final

# body functions

def run_all(sample_drive_path):
    # runs folder_mask on each subfolder in the given path
    folders = os.listdir(sample_drive_path)
    print folders
    for folder in folders:
        if folder[:3] == "cam":
            folder_mask(sample_drive_path+folder+"/", folder)

def folder_mask(input_path, output_path):
    # computes the lens smear mask for the set of pictures in a given folder
    # see slides for description of methodology
    print "Running folder: "
    print input_path
    print "Storing output in: "
    print output_path

    paths = [input_path + filename for filename in os.listdir(input_path)]
    prev_image = imageio.imread(paths[0])

    counter = 0
    for filename in paths[0:]:
        im = imageio.imread(filename)
        if counter == 0:
            # array initialization
            running_total = np.zeros(imageio.imread(filename).shape)
            subtract_total = np.zeros_like(running_total)
            gradient_running_total = np.zeros_like(running_total)
            running_average = np.zeros_like(running_total)
            total_approx_dev = np.zeros_like(running_total)
        running_total += im
        running_average = (im + running_average * counter) / (counter + 1) # recompute running average
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

# runs folder iteration function
# based on command line input

run_all(sys.argv[1])
