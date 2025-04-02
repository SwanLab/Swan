This folder contains the 7 datasets used in the project. It is also provided
the folder with the images used to build the two datasets of "Flowers". One is
the 32x32 rescaled images and the other is the dct with 3015
features truncation. 


-----Gender-----
features: 2, height and weight
groups:   2, men and women
points:   10,000

-----Microchip-----
features: 2, two test
groups:   2, rejected and accpeted
points:   118

-----4circles-----
features: 2, fictional
groups:   4, fictional
points:   1000

-----ConCircles-----
features: 2, fictional
groups:   2, fictional
points:   500

-----Iris-----
features: 4, sepal length, sepal width, petal length and petal width, 
groups:   3, setosa, versicolor and virginica
points:   500

-----MNIST-----
features: 784, the image pixels
groups:   10, hand written digits from 0 to 9
points:   10,000

-----Flowers3232-----
features: 1024, the rescaled image pixels 
groups:   5, daisy, dandelion, roses, sunflowers and tulips
points:   3451

-----FlowersDCT-----
features: 3105, feature selection from the top left corner dct matrix
groups:   5, daisy, dandelion, roses, sunflowers and tulips
points:   3451