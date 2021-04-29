# EmiyaEnginePlugin
[EmiyaEngine](https://github.com/Sg4Dylan/EmiyaEngine) 的 MATLAB/VST2 实现  

![UI](https://imgur.com/DgeYhKh.png)

### 已实现特性：
 - AkkoMode (随机抖动)  
 - CopyBand (频谱拷贝)  
 - Envelope detect (包络检测)  
 - Auto gain (自动增益调节)

### 使用方法：
 1. 准备任意支持 64 位 VST2 插件的 DAW 平台，例如 Adobe Audition CC；
 2. 自行编译或从 Release 页面下载已编译完成的插件，解压至 DAW 指定目录；
 3. 在 DAW 中找到 EmiyaEngine 效果器使用。

### 开发环境：
 - MATLAB R2021a
 - Visual Studio 2019

### Q&A：
 - Q: 使用相同的参数效果和原 Repo 不一致？
 - A: 原 Repo 中部分方法使用完整音轨遍历计算参数，在插件中以流式输入数据因此使用了近似替代；类似的，原来使用的工具实现（主要是滤波器）无法完全在 MATLAB 中找到替代，同样存在差异。

 - Q: Auto setting 开关没有反应？
 - A: 同时开启 Bypass 时程序会后台计算 Src/Dst 数值，但由于未知原因在编译的插件中数值未同步更新至界面。

