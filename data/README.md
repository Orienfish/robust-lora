# How to Prepare Dataset for LoRa Gateway Placement and End Device Assignment

 This is a tutorial walking you through the steps of how to generate dataset for our LoRa deployment problem. To prepare a dataset, you need:

* A screenshot of the area you want to deploy, e.g., `la.jpg` for our LA dataset. The latitude and longitude of the origin, which is the bottom-left corner of the image, needs to be specified in `root-of-robust-lora/main.py`.
* A file indicating the coordinates of LoRa end devices in meters with regard to the origin. For example, we use `...` for our LA dataset.
* A file indicating the coordinates of candidate gateways in meters with regard to the origin. For example, we use `...` for our LA dataset, which spreads in grids.

By the end of this tutorial, after using [GLNet](https://github.com/VITA-Group/GLNet) and the path-loss estimation algorithm proposed in [SateLoc](https://ieeexplore.ieee.org/abstract/document/9111031), you will get an path-loss matrix with the (i, j) element representing the estimated path loss in dB between the ith end device and jth gateway. Here we re-implement the path-loss estimation algorithm in SateLoc in `./path_loss_est.py`.

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

1. Copy `root-of-GLNet/entire_seg_map.npy` on Google drive into the local `root-of-robust-lora/data` folder. 

2. Run the path-loss estimation algorithm proposed in [SateLoc](https://ieeexplore.ieee.org/abstract/document/9111031) to create path loss matrix:

   ```bash
   python3 root-of-robust-lora/data/path_loss_est.py
   ```

Now you have every component to run the LoRa gateway deployment and device configuration algorithm in `root-of-robust-lora/main.py`.





