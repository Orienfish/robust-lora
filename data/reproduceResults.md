# How to reproduce LoRA experiement results 

Clone robust-lora and GLNet repositories from Github.

### Step 1: Reproduce GLNet result 

1. Open the GLNetSetup.ipynb in Google Colab.
2. Upload GLNet-master from Github to Google drive. *Make sure the first block of code in GLNetSetup.ipynb navigates to the correct location for the GLNet-master in Google drive.*
3. Modify the file "requirements.txt" as specified in GLNetSetup.ipynb.
4. Run the first 3 blocks of code.
5. Upload the test image (la.jpg) into .data/test and upload the pretrained models from the GLNet github to the saved_models folder. 
6. Run the remaining 3 blocks of code.

The result from the running the GLNet code should be uploaded to Google drive and located in the predictions folder. 

### Step 2: Convert GLNet result from png to npy file 

1. Open pngToNp.ipynb in Google Colab. *Make sure the code in block 2 navigates to where GLNet-master is located.*
2. Run the 3 blocks of code.

### Step 3: Generate PL matrix 

1. Copy "entire_seg_map.npy" created in Step 2 into robust-lora/data folder. 
2. Run the command `python3 ./alg/main.py` to generate the potential gateway and sensor locations.
3. Run the command `python3 ./data/path_loss_est.py` to create path loss matrix. 

### Step 4: Run Algorithms 

1. Run the command `python3 ./alg/main.py` again and make sure ./data/path_loss_est.py is listed as the PL matrix used in ./alg/main.py. 




