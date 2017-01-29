import imageio
import os
import numpy as np
import scipy.ndimage


def grayscale(image):
    averaged = np.mean(image, 2)
    gray = np.zeros_like(image)
    gray[:,:,0] = averaged
    gray[:,:,1] = averaged
    gray[:,:,2] = averaged
    return gray

def maxabs(values):
    return np.max(np.abs(values - values[4]))

def get_true_gradient(image):
    footprint = np.ones((3,3))
    result = np.zeros_like(image)
    result[:,:,0] = scipy.ndimage.filters.generic_filter(image[:,:,0], maxabs, footprint=footprint)
    result[:,:,1] = scipy.ndimage.filters.generic_filter(image[:,:,1], maxabs, footprint=footprint)
    result[:,:,2] = scipy.ndimage.filters.generic_filter(image[:,:,2], maxabs, footprint=footprint)
    return result

def get_prewitt_gradient(image):
    image = grayscale(image)
    footprint = np.ones((3,3))
    result = np.zeros_like(image)
    result= scipy.ndimage.prewitt(image, axis=1)
    result += scipy.ndimage.prewitt(image, axis=0)
    # result[:,:,1] = scipy.ndimage.filters.prewitt(image[:,:,1])
    # result[:,:,2] = scipy.ndimage.filters.prewitt(image[:,:,2])
    return result

def local_mean(image):
    final = np.zeros_like(image)
    final[:,:,0] = scipy.ndimage.filters.convolve(image[:,:,0], weights=np.full((5,5), 1.0/25))
    final[:,:,1] = final[:,:,0]
    final[:,:,2] = final[:,:,0]
    return final

def local_median(image):
    final = np.zeros_like(image)
    final[:,:,0] = scipy.ndimage.median_filter(image[:,:,0], (10,10))
    final[:,:,1] = final[:,:,0]
    final[:,:,2] = final[:,:,0]
    return final

def threshold(image, percent):
    threshold = (np.max(image) - np.min(image)) * percent + np.min(image)
    image[image < threshold] = 0
    image[image >= threshold] = 255
    return image

filepath = "c:/Users/Will Molter/Pictures/Geospatial vision proj1/sample_drive/cam_3/"
filename = filepath + "393408606.jpg"
running_total = np.zeros(imageio.imread(filename).shape)
subtract_total = np.zeros_like(running_total)
num_images = 40
start_image = 200

paths = [filepath + filename for filename in os.listdir(filepath)[start_image: start_image + num_images]]
prev_image = imageio.imread(paths[0])

gradient_running_total = np.zeros_like(running_total)

counter = 0
for filename in paths[0:]:
    im = imageio.imread(filename)
    running_total += im
    # gradient_running_total += get_true_gradient(im)
    print counter
    counter += 1
    subtract_total += np.abs(im - prev_image)
    prev_image = im



subtract_total_gray = grayscale(subtract_total)
subtract_normalized = threshold(subtract_total_gray, .6)
average = running_total / num_images
#
# print np.max(subtract_normalized)
# #print np.max(average)
# print np.min(subtract_normalized)

imageio.imwrite("subtract_gray.jpg", subtract_total_gray)
imageio.imwrite("subtraction.jpg", subtract_normalized)
# gradient_running_total = imageio.imread("gradient.jpg")


imageio.imwrite("average.jpg", average)
# imageio.imwrite("gradient.jpg", gradient_running_total)
# imageio.imwrite("gradient_threshold_median.jpg", local_median(threshold(gradient_running_total, .4)))