function output_value = RRMSE(image1, image2)
output_value = (sqrt((sum(sum(((abs(image1)-abs(image2)).^2))))))/(sqrt(sum(sum(abs(image1).^2))));