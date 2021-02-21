# How to Prepare Dataset for LoRa Gateway Placement and End Device Assignment

### Existing datasets

In our repo, we offer two datasets. Each dataset shares a common structure:

```
.
├── entire_seg_map.npy  // Segmented map from la.jpg using GLNet
├── gw_able.npy         // Feasibility to place gateway, 0 means the current location is
                        // over the ocean thus infeasible for gateway placement
├── gw_loc.csv          // List of candidate gateways' locations 
├── origin.csv          // Latitude and longtitude of four corners of la.jpg
├── la.jpg              // Screenshot of the satellite map of the region
├── path_loss_mat.npy   // Generated path-loss matrix
└── sr_loc.csv          // List of LoRa end devices' locations
```

More details about the two existing datasets:

* LA dataset in `./LA-dataset`: The end devices locations are generated from [PurpleAir](https://www2.purpleair.com/) air-quality sensor deployments around the LA area. The raw end devices locations downloaded from the website is in `./LA-dataset/dataLA.csv`. **(Note: origin.csv needs to be added!)**

  To generate list of LoRa end devices and candidate locations in the LA dataset, we offer a script to serve the purpose:  reading from `./LA-dataset/dataLA.csv` and generating `sr_loc.csv` and `gw_loc.csv` in the same directory:

  ```bash
  python3 ./LA-dataset/GenLocation.py
  ```

* HPWREN dataset...

### Overview

This tutorial will you through the steps of how to generate dataset for our LoRa deployment problem.

To prepare a dataset, you need:

* A screenshot of the area you want to deploy, e.g., `./LA-dataset/la.jpg` for our LA dataset. The latitude and longitude of the origin is also needed, which is stored in a file.
* A file indicating the coordinates of LoRa end devices in meters with regard to the origin. For example, we use `./LA-dataset/sr_loc.csv` for our LA dataset.
* A file indicating the coordinates of candidate gateways in meters with regard to the origin. For example, we use `./LA-dataset/gw_loc.csv` for our LA datase.

By the end of this tutorial, using [GLNet](https://github.com/VITA-Group/GLNet) and the path-loss estimation algorithm proposed in [SateLoc](https://ieeexplore.ieee.org/abstract/document/9111031), you will get a path-loss matrix (e.g. `./LA-dataset/path_loss_mat.npy`) . The (i, j)th element of the matrix represents the estimated path loss in dB between the ith end device and jth gateway. Here we re-implement the path-loss estimation algorithm in SateLoc in `./path_loss_est.py`.

### Step 0: Preparation 

There are two options to reproduce GLNet result: running locally or on cloud. In our case, we use Google Colab to avoid complicated settings of GPU. We present the instructions of running GLNet in Google colab in the following lines:

1. Clone robust-lora (this repo) and [GLNet repo](https://github.com/VITA-Group/GLNet) from Github. 

   ```bash
   git clone https://github.com/Orienfish/robust-lora.git
   git clone https://github.com/VITA-Group/GLNet.git
   ```

2. Navigate to the root directory of GLNet and modify `requirements.txt` to:

```
numpy
torch
torchvision
tqdm
tensorboardX
Pillow==6.2.2
opencv-python
```

3. Download the pretrained models from `GLNet` and put them into folder `root-of-GLNet/saved_models`. This step is the same as the Evaluation section in the [GLNet tutorial](https://github.com/VITA-Group/GLNet).
   * [fpn_deepglobe_global.pth](https://drive.google.com/file/d/1xUJoNEzj5LeclH9tHXZ2VsEI9LpC77kQ/view?usp=sharing)
   * [fpn_deepglobe_global2local.pth](https://drive.google.com/file/d/1_lCzi2KIygcrRcvBJ31G3cBwAMibn_AS/view?usp=sharing)
   * [fpn_deepglobe_local2global.pth](https://drive.google.com/file/d/198EcAO7VN8Ujn4N4FBg3sRgb8R_UKhYv/view?usp=sharing)
4. Upload the whole directory of `GLNet` to your Google drive. 
5. Upload the image you want to test into Google drive folder `root-of-GLNet/data/test`, for example, `la.jpg` for the LA dataset in our experiments.

### Step 1: Reproduce GLNet results and store in `.npy` file

1. Open the `./GLNetSetup.ipynb` in Google Colab. **Make sure the first block of code in `GLNetSetup.ipynb` navigates to the correct location for the`GLNet` in Google drive.** In our case, we upload `GLNet` to the root of the drive, thus the following lines will do the job of mounting drive and installing required packages:

   ```python
   drive.mount('/content/drive')
   %cd drive/MyDrive/GLNet
   %pip install -r requirements.txt
   ```

2. Run the blocks in `./GLNetSetup.ipynb` one by one. You need to specify the name of the image you want to test in the 2nd line of the 3rd block.

The segmented images resulted from running GLNet should be uploaded to Google drive and located in the `root-of-GLNet/prediction` folder. In our case, it is named `root-of-GLNet/prediction/sat_test_entire_mask.png`. 

In the last block in `./GLNetSetup.ipynb` converts the generated image to an array and stores in `root-of-GLNet/entire_seg_map.npy`.

### Step 2: Generate path-loss matrix 

1. Download `root-of-GLNet/entire_seg_map.npy` on Google drive to the local `root-of-robust-lora/data` folder. 

2. Run the path-loss estimation algorithm proposed in [SateLoc](https://ieeexplore.ieee.org/abstract/document/9111031) to create path loss matrix. Make sure you have specified the paths to (i) device locations, (ii) candidate gateway locations and (iii) segmented maps in `path_loss_est.py`. **(Note: it is recommended to add command line arguments to specify these paths)**

   ```bash
   python3 root-of-robust-lora/data/path_loss_est.py
   ```

Now you have every component to run the LoRa gateway deployment and device configuration algorithm in `root-of-robust-lora/alg/main.py`. Make sure you have specified the path to (i) device locations, (ii) candidate gateway locations and (iii) path-loss matrix files in class `DataParams` (line 74) in `main.py`.
