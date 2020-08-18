# robust-lora

Use the [lorawan](https://github.com/signetlabdei/lorawan/commit/c2833b75662e7cdadd78b10a2a1d6d8f3c6a5769) module developed by D. Margrin et al.

**Changes to be made in our simulation**:

* To enable confirmed traffic, change line 343 of `lorawan/model/end-device-lorawan-mac.h` to as follows:

  ```c++
  bool waitingAck = true;
  ```

* 

