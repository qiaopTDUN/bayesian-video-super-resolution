		`Bayesian video super-resolution`
		
		This is Matlab implementation of a Bayesian video super-resolution method [1] based on [2]. 
Compared with the implementation in [2], this implementation is more closer to the details 
described in [1], e.g., the update of optical flow.
		
		`Description of folders`
		
		/data: contains test sequneces, and the results produced by [1]. The original data can be 
downloaded from [3].
		/celiu_optical_flow: contains optical flow estimation codes, which be downloaded from [4].
		
		`Usage`
		1. Before running, please check and change the path for celiu_optical_flow!
		2. run mfsr_cvpr2011_main.m
			2.1 city sequence as default, recover the 16th frame in upscale=2 scenario.
			2.2 The recovered image is stored in variable TP.Ik
		3. You can change the parameters setting in initParam.m, and run mfsr_cvpr2011_main.m. 
Different parameters result in different image recover quality. Just enjoy parameters tunning, and have a fun!

		`Contact`
		Peng Qiao, Email: pengqiao@nudt.edu.cn
		
		`Citation`
		[1] C. Liu, and D. Sun, "On Bayesian Adaptive Video Super Resolution," IEEE Trans. on Pattern Analysis and Machine Intelligence (PAMI), Feb. 2014. 
		[2] https://github.com/seunghwanyoo/bayesian_vid_sr
		[3] http://people.csail.mit.edu/celiu/CVPR2011/default.html
		[4] http://people.csail.mit.edu/celiu/OpticalFlow/
