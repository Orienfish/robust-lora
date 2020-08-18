# robust-lora

We use the [lorawan](https://github.com/Orienfish/lorawan) module forked from the [original implementation](https://github.com/signetlabdei/lorawan) by D. Margrin et al.

**Following changes are made in our simulation**:

* To enable confirmed traffic, change line 343 of `lorawan/model/end-device-lorawan-mac.h` to as follows:

  ```c++
  bool waitingAck = true;
  ```

* 

