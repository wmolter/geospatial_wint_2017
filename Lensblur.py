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

def local_median(image, size=5):
    final = np.zeros_like(image)
    # final[:,:,0] = scipy.ndimage.median_filter(image[:,:,0], (10,10))
    # final[:,:,1] = final[:,:,0]
    # final[:,:,2] = final[:,:,0]
    return each_layer_same(scipy.ndimage.median_filter(image[:,:,0], (size,size)), final)

def threshold(image, percent):
    threshold = (np.max(image) - np.min(image)) * percent + np.min(image)
    image[image < threshold] = 0
    image[image >= threshold] = 255
    return image

# filepath depends on who's running the code

filepath_will = "c:/Users/Will Molter/Pictures/Geospatial vision proj1/sample_drive/cam_3/"
filepath_jacob = "/Users/jdbruce/Downloads/WQ2017/Geospatial/sample_drive/cam_3/"

filepath = filepath_jacob # change this line when you change users


# loop portion
# comment here so it doesn't loop


# num_images = 1000
# start_image = 100

# paths = [filepath + filename for filename in os.listdir(filepath)[start_image: start_image + num_images]]
# prev_image = imageio.imread(paths[0])

# counter = 0
# for filename in paths[0:]:
#     im = imageio.imread(filename)
#     if counter == 0:
#         running_total = np.zeros(imageio.imread(filename).shape)
#         subtract_total = np.zeros_like(running_total)
#         gradient_running_total = np.zeros_like(running_total)
#     running_total += im
#     # gradient_running_total += get_true_gradient(im)
#     print counter
#     counter += 1
#     subtract_total += np.abs(im - prev_image)
#     prev_image = im



# subtract_total_gray = grayscale(subtract_total)
# subtract_normalized = threshold(subtract_total_gray, .6)
# average = running_total / num_images
# print np.min(subtract_normalized)

# imageio.imwrite("subtract_gray.jpg", subtract_total_gray)
# imageio.imwrite("subtraction.jpg", subtract_normalized)
# gradient_running_total = imageio.imread("gradient.jpg")

average = imageio.imread("average.jpg")
# average_sobel = sobel_gradient(average)
average_prewitt = imageio.imread("average_prewitt.jpg")
# average_grad = get_true_gradient(average)
# average_prewitt = get_prewitt_gradient(average)
average_sobel = imageio.imread("average_sobel.jpg")

average_prewitt_mean = local_mean(average_prewitt, 10)
# average_prewitt_mean_mean = local_mean(average_prewitt_mean)
average_prewitt_median_5 = local_median(average_prewitt, 5)
average_prewitt_median_10 = local_median(average_prewitt, 10)
average_prewitt_median_15 = local_median(average_prewitt, 15)
average_prewitt_median_20 = local_median(average_prewitt, 20)

imageio.imwrite("average.jpg", average)
imageio.imwrite("average_prewitt.jpg", average_prewitt)
imageio.imwrite("average_sobel.jpg", average_sobel)
imageio.imwrite("average_prewitt_mean.jpg", average_prewitt_mean)
imageio.imwrite("average_prewitt_median_5.jpg", average_prewitt_median_5)
imageio.imwrite("average_prewitt_median_10.jpg", average_prewitt_median_10)
imageio.imwrite("average_prewitt_median_15.jpg", average_prewitt_median_15)
imageio.imwrite("average_prewitt_median_20.jpg", average_prewitt_median_20)
# imageio.imwrite("average_prewitt_mean_mean.jpg", average_prewitt_mean_mean)

# imageio.imwrite("gradient.jpg", gradient_running_total)
# imageio.imwrite("gradient_threshold_median.jpg", local_median(threshold(gradient_running_total, .4)))






