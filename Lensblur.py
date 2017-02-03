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
filepath_jacob = "/Users/jdbruce/Downloads/WQ2017/Geospatial/sample_drive/cam_0/"

filepath = filepath_jacob # change this line when you change users


# loop portion
# comment here so it doesn't loop

    def folder_mask(path):


    num_images = 10
    start_image = 100

    paths = [filepath + filename for filename in os.listdir(filepath)[start_image: start_image + num_images]]
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

    # get mask 1 (from grad)

    grad_norm = normalize(gradient_running_total)
    imageio.imwrite("grad_norm.jpg", grad_norm)
    grad_norm_threshold = threshold(grad_norm, 40)
    imageio.imwrite("grad_norm_threshold.jpg", grad_norm_threshold)
    grad_norm_threshold_median = local_median(grad_norm_threshold)
    imageio.imwrite("grad_norm_threshold_median.jpg", grad_norm_threshold_median)

    mask_1_final = grad_norm_threshold_median

    # get mask 2 (from dev)

    dev_normalized = normalize(total_approx_dev)
    imageio.imwrite("deviation.jpg", dev_normalized)

    # dev_normalized = imageio.imread("deviation.jpg")
    # imageio.imwrite("dev_normalized.jpg", dev_normalized)
    # # dev_normalized_mean = local_mean(dev_normalized, 50)
    dev_normalized_median =  local_median(dev_normalized)  #imageio.imread("dev_normalized_median.jpg") #
    # dev_normalized_median_median = local_median(dev_normalized_median, 20)  #
    # imageio.imwrite("dev_normalized_median_median.jpg", dev_normalized_median_median)
    imageio.imwrite("dev_normalized_median.jpg", dev_normalized_median)

    # dev_normalized_median_prewitt = prewitt_gradient(dev_normalized_median)
    # imageio.imwrite("dev_normalized_median_prewitt.jpg", dev_normalized_median_prewitt)


    # dev_normalized_median_sharpen = sharpen(dev_normalized_median)
    # imageio.imwrite("dev_normalized_median_sharpen.jpg", dev_normalized_median_sharpen)


    #
    dev_normalized_median_threshold = threshold(dev_normalized_median, 10)
    imageio.imwrite("dev_normalized_median_threshold.jpg", dev_normalized_median_threshold)

    mask_2_final = dev_normalized_median_threshold

    # final output_array

    # mask_1_final = imageio.imread("grad_norm_threshold_median.jpg")
    # mask_2_final = imageio.imread("dev_normalized_median_threshold.jpg")

    final_mask = NAND(mask_1_final, mask_2_final)
    imageio.imwrite("final_mask.jpg", final_mask)





    # dev_normalized_median_gradient = true_gradient(dev_normalized_median)
    # imageio.imwrite("dev_normalized_median_gradient.jpg", dev_normalized_median_gradient)

    # dev_normalized_prewitt = prewitt_gradient(dev_normalized)
    # imageio.imwrite("dev_normalized_prewitt.jpg",dev_normalized_prewitt)

    # average = imageio.imread("average.jpg")
    # # average_sobel = sobel_gradient(average)
    # average_prewitt = imageio.imread("average_prewitt.jpg")
    # # average_grad = get_true_gradient(average)
    # # average_prewitt = get_prewitt_gradient(average)
    # average_sobel = imageio.imread("average_sobel.jpg")
    #
    # average_prewitt_mean = local_mean(average_prewitt, 10)
    # # average_prewitt_mean_mean = local_mean(average_prewitt_mean)
    # average_prewitt_median_5 = local_median(average_prewitt, 5)
    # average_prewitt_median_10 = local_median(average_prewitt, 10)
    # average_prewitt_median_15 = local_median(average_prewitt, 15)
    # average_prewitt_median_20 = local_median(average_prewitt, 20)
    #
    # imageio.imwrite("average.jpg", average)
    # imageio.imwrite("average_prewitt.jpg", average_prewitt)
    # imageio.imwrite("average_sobel.jpg", average_sobel)
    # imageio.imwrite("average_prewitt_mean.jpg", average_prewitt_mean)
    # imageio.imwrite("average_prewitt_median_5.jpg", average_prewitt_median_5)
    # imageio.imwrite("average_prewitt_median_10.jpg", average_prewitt_median_10)
    # imageio.imwrite("average_prewitt_median_15.jpg", average_prewitt_median_15)
    # imageio.imwrite("average_prewitt_median_20.jpg", average_prewitt_median_20)
    # imageio.imwrite("average_prewitt_mean_mean.jpg", average_prewitt_mean_mean)

    # imageio.imwrite("gradient.jpg", gradient_running_total)
