
# 前言
本文为笔者项目前处理所需算法，仅供参考学习，请遵顼GPL-3.0协议，转载请备注署名。
> This project is licensed under the GNU General Public License v3.0 (GPL-3.0).  See the [LICENSE ](https://www.gnu.org/licenses/gpl-3.0.txt)file for details.

---
FPGA图像处理中常用的图像缩放算法，包括最近邻插值算法、双线性邻插值算法、双三次邻插值算法，以及双线性的MATLAB与FPGA实现。

# 双线性邻域插值算法

## 基础理论
**核心思想：** 标准的线性插值是在一维直线上，根据已知两点的坐标和值，来计算中间某个未知点的值。而双线性插值指的是做了两次线性插值，是在一个二维平面上的。而 **双线性插值 (Bilinear Interpolation)**，顾名思义，就是做了两次线性插值，但它是在一个二维平面上。当我们要计算目标图像中某个点 `(x, y)` 的像素值时，我们先找到它对应在源图像中的位置 `(src_x, src_y)`。这个位置通常是浮点数，落在四个真实像素点之间。
![image.png](https://obsidian-picturebed-1256135654.cos.ap-nanjing.myqcloud.com/obsidion/20250725153727827.png?imageSlim)
**算法步骤：**
1. **坐标映射**：将目标图像的坐标 `(dst_x, dst_y)` 映射回源图像的坐标 `(src_x, src_y)`。
2. **找到四个邻近点**：根据 `(src_x, src_y)` 找到它周围的四个像素点，即上图的 `Q11, Q12, Q21, Q22`。（如上图）
3. **第一次线性插值**：在 x 方向上，分别对上下两行进行线性插值，计算出两个临时点 `R1` 和 `R2` 的值。
4. **第二次线性插值**：在 y 方向上，对 `R1` 和 `R2` 进行线性插值，最终得到 `P` 点的像素值。

**公式：** 


$$
f(P) = (1 - d_x)(1 - d_y) f(Q_{11}) + d_x(1 - d_y) f(Q_{21}) + (1 - d_x)d_y f(Q_{12}) + d_x d_y f(Q_{22})
$$

**原理实例化：**
![[ISP之常用图像缩放算法 2025-07-25 20.53.47.excalidraw|1000]]

## 笔者理论改进
笔者最终目标是实现再FPGA可以高效运行的双线性插值，为此需要在一开始从算法层面去改进。出于习惯，个人觉得以下公式更清晰：
$$
I(i, j) = [1 - y(i)] \times \left\{ [(1 - x(j)) \times T_A(r) + x(j) \times T_B(r)] \right\} + y(i) \times \left\{ [(1 - x(j)) \times T_C(r) + x(j) \times T_D(r)] \right\}
$$
为了清晰地表达，约定以下符号：
- **Ws​,Hs​**：源图像的宽度和高度。
- **Wd​,Hd​**：目标图像的宽度和高度。
- **(dstx​,dsty​)**：目标图像中的像素坐标。
- **(srcx​,srcy​)**：通过映射在源图像中得到的浮点坐标。
- **(x1​,y1​)**：包围 (srcx​,srcy​) 的四点矩阵的左上角整数坐标，即 x1​=⌊srcx​⌋,y1​=⌊srcy​⌋。
- **dx,dy**：浮点坐标的小数部分，即 dx=srcx​−x1​,dy=srcy​−y1​。
- **N**：定点运算的小数位数，这里 N=10。
- **S**：缩放因子 (Scale)，S=2N=1024。
- **Xfixed​**：表示变量 X 的定点整数形式。
- **TA,TB,TC,TD**：分别代表源图像中四个邻近像素 f(x1​,y1​), f(x1​+1,y1​), f(x1​,y1​+1), f(x1​+1,y1​+1) 的值。
- **≪N**：代表逻辑左移N位（乘以 S=2N）。
- **≫N**：代表算术右移N位（除以 S=2N 并向下取整）。


**改进1：**
小数在FPGA中是难以接受的，所以将1扩大1024倍，为此得到
$$
I(i, j) = [1024 - y(i)] \times \left\{ [(1024 - x(j)) \times T_A(r) + x(j) \times T_B(r)] \right\} + y(i) \times \left\{ [(1024 - x(j)) \times T_C(r) + x(j) \times T_D(r)] \right\}
$$

**改进2：**
上述公式存在6次乘法，而且在实际中位宽较大，比较消耗资源，为此采取差分方式来代替上述表达式，差分情况下只需要做三次乘法，对于1024的除法采用移位即可。
$$
I(i,j)\;=\;
\bigl(1024-y(i)\bigr)\Bigl[T_A(r)+\bigl(T_B(r)-T_A(r)\bigr)\tfrac{x(j)}{1024}\Bigr]
\;+\;
y(i)\Bigl[T_C(r)+\bigl(T_D(r)-T_C(r)\bigr)\tfrac{x(j)}{1024}\Bigr]
$$
## 分步数学推导
### 前处理 / 地址生成单元
这部分的目标是根据目标像素坐标 `(dst_x, dst_y)`，计算出插值所需的四个源像素位置和两个权重。
**1. 坐标映射：**
将目标坐标范围 `[0, Wd-1]` 映射到源坐标范围 `[0, Ws-1]`。
标准浮点公式：
$$
\begin{align*}
src_x &= dst_x \cdot \frac{W_s - 1}{W_d - 1} \\
src_y &= dst_y \cdot \frac{H_s - 1}{H_d - 1}
\end{align*}
$$
定点实现公式:
$$
\begin{align*}
src\_x_{fixed} &= \left\lfloor \frac{dst_x \cdot (W_s - 1) \cdot S}{W_d - 1} \right\rfloor \\
src\_y_{fixed} &= \left\lfloor \frac{dst_y \cdot (H_s - 1) \cdot S}{H_d - 1} \right\rfloor
\end{align*}
$$
**2. 分解整数坐标与小数权重**
从定点坐标 `src_x_fixed` 中分离出整数部分 `x1` 和归一化的小数权重 `x_norm`。
整数坐标 `x1,y1` 公式： 
$$
\begin{align*}
x_1 &= \left\lfloor \frac{src\_x_{fixed}}{S} \right\rfloor \quad \Leftrightarrow \quad x_1 = src\_x_{fixed} \gg N \\
y_1 &= \left\lfloor \frac{src\_y_{fixed}}{S} \right\rfloor \quad \Leftrightarrow \quad y_1 = src\_y_{fixed} \gg N
\end{align*}
$$
小数坐标`x_norm` ,`y_norm`公式：
$$
\begin{align*}
dx &= src_x - x_1 \implies dx_{fixed} = (src_x - x_1) \cdot S = src\_x_{fixed} - x_1 \cdot S \\
x_{norm} &= src\_x_{fixed} \pmod{S} \quad \Leftrightarrow \quad x_{norm} = src\_x_{fixed} \,\&\, (S - 1)
\end{align*}
$$

`y_norm` 同理
### 插值计算核心
省流：为什么这么做？因为这个版本逻辑对于fpga来说既能保持精度（与浮点存在±1误差），又符合计算逻辑
```
// R1 = TA << 10 + x * (TB - TA)
// R2 = TC << 10 + x * (TD - TC)
// I  = R1 >> 10 + [y * (R2 - R1)] >> 20
```
**1. X方向的线性插值**
计算中间点 R1​ (在 y1​ 行) 和 R2​ (在 y1​+1 行) 的值。
标准浮点公式：
$$
\begin{align*}
R_1 &= (1 - dx) \cdot TA + dx \cdot TB \\
R_2 &= (1 - dx) \cdot TC + dx \cdot TD
\end{align*}
$$
定点公式：
$$
R_{1_{fixed}} = (S - x_{norm}) \cdot TA + x_{norm} \cdot TB
$$
展开
$$ 
R\_1_{fixed} = S \cdot TA - x_{norm} \cdot TA + x_{norm} \cdot TB = (TA \ll N) + x_{norm} \cdot (TB - TA)
$$
R_2fix同理
**2. Y方向的线性插值**
利用 R1​ 和 R2​ 计算最终像素值 I。
标准浮点公式：
$$
I = (1 - dy) \cdot R_1 + dy \cdot R_2
$$
定点实现公式:
$$
I_{final} = \left\lfloor \frac{(S - y_{norm}) \cdot R_{1_{fixed}} + y_{norm} \cdot R_{2_{fixed}}}{S^2} \right\rfloor
$$
拆分展开：
$$
\begin{align*}
I_{final} &= \left\lfloor \frac{S \cdot R_{1_{fixed}} + y_{norm} \cdot (R_{2_{fixed}} - R_{1_{fixed}})}{S^2} \right\rfloor \\
I_{final} &= \left\lfloor \frac{R_{1_{fixed}}}{S} + \frac{y_{norm} \cdot (R_{2_{fixed}} - R_{1_{fixed}})}{S^2} \right\rfloor
\end{align*}
$$
**3. 最终步骤**
$$
I_{final} = (R_{1_{fixed}} \gg N) + \left( (y_{norm} \cdot (R_{2_{fixed}} - R_{1_{fixed}})) \gg 2N \right)
$$
## python定点实现
### code
```python
import cv2  
import numpy as np  
import numpy as ny  
import matplotlib.pyplot as plt  
from scipy.stats import yeojohnson_normmax  
  
# --- Matplotlib 中文显示设置 ---import matplotlib  
  
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体  
matplotlib.rcParams['axes.unicode_minus'] = False  # 正常显示负号  
  
def bilinear_verilog(src_img,dst_shape):  
    FRAC_BITS = 10  
    SCALE = 1 << FRAC_BITS # 1024  
  
    src_h,src_w = src_img.shape[:2]  
    dst_h,dst_w = dst_shape  
  
    if len(src_img.shape) == 3:  
        channels = src_img.shape[2]  
    else:  
        channels =1  
  
    dst_img = np.zeros((dst_h,dst_w,channels),np.uint8)  
  
    for dst_y in range(dst_h):  
        for dst_x in range(dst_w):  
            # 坐标映射  
            src_x_fixed = (dst_x * (src_w-1) * SCALE)//(dst_w-1)  
            src_h_fixed = (dst_y * (src_h-1) * SCALE)//(dst_h-1)  
  
            # 分解坐标和权重  
            x1 = src_x_fixed >> FRAC_BITS  
            y1 = src_h_fixed >> FRAC_BITS  
            x_norm = src_x_fixed & (SCALE-1)  
            y_norm = src_h_fixed & (SCALE-1)  
  
            # 边界处理  
            x2 = x1 + 1  
            y2 = y1 + 1  
            x1 = max(0, x1);x2 = min(src_w - 1, x2)  
            y1 = max(0, y1);y2 = min(src_h - 1, y2)  
  
            # 核心公式  
            # R1 = TA << 10 + x * (TB - TA)  
            # R2 = TC << 10 + x * (TD - TC)            # I = R1 >> 10 + [y * (R2 - R1)] >> 20  
            for c in range(channels):  
                # 获取4个临近点  
                TA = int(src_img[y1,x1,c])  
                TB = int(src_img[y1,x2,c])  
                TC = int(src_img[y2, x1, c])  
                TD = int(src_img[y2, x2, c])  
  
                # x方向插值  
                diff1 = TB - TA  
                diff2 = TD - TC  
                R1 = (TA << FRAC_BITS) + x_norm * diff1  
                R2 = (TB << FRAC_BITS) + x_norm * diff2  
  
                # y方向插值  
                term1 = R1 >> FRAC_BITS  
  
                diff_R = R2 - R1  
                mult_v = y_norm * diff_R  
  
                rounding_offset = 1 << (2 * FRAC_BITS - 1)  
                term2 = (mult_v + rounding_offset) >> (2 * FRAC_BITS)  
  
                final_pixel_val = term1 + term2  
  
                dst_img[dst_y, dst_x, c] = np.clip(final_pixel_val,0,255)  
  
        if channels == 1:  
            dst_img = np.squeeze(dst_img, axis=2)  
  
    return dst_img  
  
# mian  
if __name__=='__main__':  
    img = cv2.imread('tenshu.jpg')  
    img = cv2.cvtColor(img,cv2.COLOR_BGR2RGB)  
  
  
    target_shape = (480,640) # 高度、宽度  
  
    zoomed_verilog = bilinear_verilog(img,target_shape)  
  
    target_size_cv2 = (target_shape[1], target_shape[0])  
    zoomed_cv2 = cv2.resize(img,target_size_cv2,interpolation=cv2.INTER_LINEAR)  
  
    # 显示结果  
    plt.figure(figsize=(15, 7))  
    plt.suptitle("Verilog风格定点实现 vs OpenCV实现", fontsize=16)  
  
    plt.subplot(1, 3, 1)  
    plt.imshow(img)  
    plt.title(f"原始图像 ({img.shape[1]}x{img.shape[0]})")  
  
    plt.subplot(1, 3, 2)  
    plt.imshow(zoomed_verilog)  
    plt.title(f"Verilog风格定点实现 ({target_shape[1]}x{target_shape[0]})")  
  
    plt.subplot(1, 3, 3)  
    plt.imshow(zoomed_cv2)  
    plt.title(f"OpenCV cv2.resize")  
  
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  
    plt.show()
```


### result
![image.png](https://obsidian-picturebed-1256135654.cos.ap-nanjing.myqcloud.com/obsidion/20250726191459555.png?imageSlim)

### 优化
除法计算过于多，每个像素都进行一次除法，这是极其反优化的。为此优化只做一次，预先计算一个高精度的“组合缩放因子”，这个因子包含了原始的缩放比例和用于最终结果的SCALE。
```python
# -*- coding: utf-8 -*-  
import cv2  
import numpy as np  
import matplotlib.pyplot as plt  
  
# --- Matplotlib 中文显示设置 ---import matplotlib  
matplotlib.rcParams['font.sans-serif'] = ['SimHei']  # 使用黑体  
matplotlib.rcParams['axes.unicode_minus'] = False  # 正常显示负号  
  
def bilinear_verilog_optimized(src_img, dst_shape):  
    # --- 定点数和精度定义 ---    FRAC_BITS = 10  
    SCALE = 1 << FRAC_BITS  # 1024  
  
    # --- 图像基本信息 ---    src_h, src_w = src_img.shape[:2]  
    dst_h, dst_w = dst_shape  

    if len(src_img.shape) == 3:  
        channels = src_img.shape[2]  
    else:  
        channels = 1  
        # 统一将灰度图也处理为3维，方便计算  
        src_img = np.expand_dims(src_img, axis=-1)  
  
    dst_img = np.zeros((dst_h, dst_w, channels), np.uint8)  
  

    RATIO_PRECISION_BITS = 16  
    RATIO_SCALE = 1 << RATIO_PRECISION_BITS  
  
    if dst_w > 1:  
        combined_scale_x = ((src_w - 1) * SCALE * RATIO_SCALE) // (dst_w - 1)  
    else:  
        combined_scale_x = 0  
  
    if dst_h > 1:  
        combined_scale_y = ((src_h - 1) * SCALE * RATIO_SCALE) // (dst_h - 1)  
    else:  
        combined_scale_y = 0  
  
    # --- 遍历目标图像 ---    for dst_y in range(dst_h):  
        for dst_x in range(dst_w):  
            # --- 坐标映射 (优化后) ---  
            src_x_fixed = (dst_x * combined_scale_x) >> RATIO_PRECISION_BITS  
            src_y_fixed = (dst_y * combined_scale_y) >> RATIO_PRECISION_BITS  
  
            # --- 分解坐标和权重 ---            x1 = src_x_fixed >> FRAC_BITS  
            y1 = src_y_fixed >> FRAC_BITS  
            x_norm = src_x_fixed & (SCALE - 1)  
            y_norm = src_y_fixed & (SCALE - 1)  
  
            # --- 边界处理 ---            x2 = x1 + 1  
            y2 = y1 + 1  
            x1 = max(0, x1); x2 = min(src_w - 1, x2)  
            y1 = max(0, y1); y2 = min(src_h - 1, y2)  
  
            # --- 核心插值公式 ---            # R1 = (TA << 10) + x_norm * (TB - TA)            # R2 = (TC << 10) + x_norm * (TD - TC)            # I  = (R1 >> 10) + [y_norm * (R2 - R1)] >> 20            for c in range(channels):  
                # 获取4个临近点  
                TA = int(src_img[y1, x1, c])  
                TB = int(src_img[y1, x2, c])  
                TC = int(src_img[y2, x1, c])  
                TD = int(src_img[y2, x2, c])  
  
                # X方向插值  
                diff1 = TB - TA  
                diff2 = TD - TC  
                R1 = (TA << FRAC_BITS) + x_norm * diff1  
                # *** 修正点: R2的计算应基于TC和TD ***  
                R2 = (TC << FRAC_BITS) + x_norm * diff2  
  
                # Y方向插值  
                term1 = R1 >> FRAC_BITS  
                diff_R = R2 - R1  
                mult_v = y_norm * diff_R  
  
                # 添加舍入值以提高精度  
                rounding_offset = 1 << (2 * FRAC_BITS - 1)  
                term2 = (mult_v + rounding_offset) >> (2 * FRAC_BITS)  
  
                final_pixel_val = term1 + term2  
  
                dst_img[dst_y, dst_x, c] = np.clip(final_pixel_val, 0, 255)  
  
    if channels == 1:  
        dst_img = np.squeeze(dst_img, axis=-1)  
  
    return dst_img  
  
# --- 主程序 ---if __name__ == '__main__':  
    try:  
        img = cv2.imread('tenshu.jpg')  
        img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)  
    except Exception as e:  
        print(f"无法加载图像 'tenshu.jpg'，请确保文件存在。错误: {e}")  
        # 创建一个备用测试图像  
        img = np.zeros((100, 100, 3), dtype=np.uint8)  
        img[0:50, 0:50] = [255, 0, 0]  
        img[50:100, 0:50] = [0, 255, 0]  
        img[0:50, 50:100] = [0, 0, 255]  
        img[50:100, 50:100] = [255, 255, 0]  
  
  
    target_shape = (480, 640)  # 高度、宽度  
  
    zoomed_verilog = bilinear_verilog_optimized(img, target_shape)  
  
    target_size_cv2 = (target_shape[1], target_shape[0])  
    zoomed_cv2 = cv2.resize(img, target_size_cv2, interpolation=cv2.INTER_LINEAR)  
  
    # 显示结果  
    plt.figure(figsize=(15, 7))  
    plt.suptitle("优化后的定点实现 vs OpenCV实现", fontsize=16)  
  
    plt.subplot(1, 3, 1)  
    plt.imshow(img)  
    plt.title(f"原始图像 ({img.shape[1]}x{img.shape[0]})")  
  
    plt.subplot(1, 3, 2)  
    plt.imshow(zoomed_verilog)  
    plt.title(f"优化后的定点实现 ({target_shape[1]}x{target_shape[0]})")  
  
    plt.subplot(1, 3, 3)  
    plt.imshow(zoomed_cv2)  
    plt.title(f"OpenCV cv2.resize")  
  
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])  
    plt.show()
```


## FPGA实现
### 模块设计
![[ISP之常用图像缩放算法 2025-07-26 19.32.04.excalidraw|2000]]


#### Verilog实现
参考：韩冰 《基于MATLAB和FPGA的图像处理教程》
（自己写一周，处理快慢问题搞不定。。。）
```Verilog
module bilinear_interpolation
#(
    // 输入改为 1024×768，输出改为 640×480
        parameter          C_SRC_IMG_WIDTH  = 11'd1024             ,
        parameter          C_SRC_IMG_HEIGHT = 11'd768              ,
        parameter          C_DST_IMG_WIDTH  = 11'd640              ,
        parameter          C_DST_IMG_HEIGHT = 11'd480              ,
    // ratio = floor(1024/640 * 2^16) = floor(1.6 * 65536) = 104857
        parameter          C_X_RATIO        = 18'd104857           ,  // 扩宽到 18 位以免溢出
        parameter          C_Y_RATIO        = 18'd104857           
)
(
        input  wire                                       clk_in1             ,
        input  wire                                       clk_in2             ,
        input  wire                                       rst_n               ,
    
    //  Image data prepared to be processed
        input  wire                                       per_img_vsync       ,//  Prepared Image data vsync valid signal
        input  wire                                       per_img_href        ,//  Prepared Image data href vaild  signal
        input  wire         [               7: 0]         per_img_gray        ,//  Prepared Image brightness input
    
    //  Image data has been processed
        output reg                                        post_img_vsync      ,//  processed Image data vsync valid signal
        output reg                                        post_img_href       ,//  processed Image data href vaild  signal
        output reg          [               7: 0]         post_img_gray        //  processed Image brightness output
);
//----------------------------------------------------------------------
reg                                      per_img_href_dly       ;

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        per_img_href_dly <= 1'b0;
    else
        per_img_href_dly <= per_img_href;
end

wire                                     per_img_href_neg       ;

assign   per_img_href_neg = per_img_href_dly & ~per_img_href;

reg            [              10: 0]     img_vs_cnt             ;//  from 0 to C_SRC_IMG_HEIGHT - 1

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        img_vs_cnt <= 11'b0;
    else
    begin
        if(per_img_vsync == 1'b0)
            img_vs_cnt <= 11'b0;
        else
        begin
            if(per_img_href_neg == 1'b1)
                img_vs_cnt <= img_vs_cnt + 1'b1;
            else
                img_vs_cnt <= img_vs_cnt;
        end
    end
end

reg            [              10: 0]     img_hs_cnt             ;//  from 0 to C_SRC_IMG_WIDTH - 1

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        img_hs_cnt <= 11'b0;
    else
    begin
        if((per_img_vsync == 1'b1)&&(per_img_href == 1'b1))
            img_hs_cnt <= img_hs_cnt + 1'b1;
        else
            img_hs_cnt <= 11'b0;
    end
end

//----------------------------------------------------------------------
reg            [               7: 0]     bram_a_wdata           ;

always @(posedge clk_in1)
begin
    bram_a_wdata <= per_img_gray;
end

reg            [              11: 0]     bram_a_waddr           ;

always @(posedge clk_in1)
begin
    bram_a_waddr <= {img_vs_cnt[2:1],10'b0} + img_hs_cnt;
end

reg                                      bram1_a_wenb           ;

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        bram1_a_wenb <= 1'b0;
    else
        bram1_a_wenb <= per_img_vsync & per_img_href & ~img_vs_cnt[0];
end

reg                                      bram2_a_wenb           ;

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        bram2_a_wenb <= 1'b0;
    else
        bram2_a_wenb <= per_img_vsync & per_img_href & img_vs_cnt[0];
end

reg            [              10: 0]     fifo_wdata             ;

always @(posedge clk_in1)
begin
    fifo_wdata <= img_vs_cnt;
end

reg                                      fifo_wenb              ;

always @(posedge clk_in1)
begin
    if(rst_n == 1'b0)
        fifo_wenb <= 1'b0;
    else
    begin
        if((per_img_vsync == 1'b1)&&(per_img_href == 1'b1)&&(img_hs_cnt == C_SRC_IMG_WIDTH - 1'b1))
            fifo_wenb <= 1'b1;
        else
            fifo_wenb <= 1'b0;
    end
end

//----------------------------------------------------------------------
//  bram & fifo rw
reg            [              11: 0]     even_bram1_b_raddr     ;
reg            [              11: 0]     odd_bram1_b_raddr      ;
reg            [              11: 0]     even_bram2_b_raddr     ;
reg            [              11: 0]     odd_bram2_b_raddr      ;
wire           [               7: 0]     even_bram1_b_rdata     ;
wire           [               7: 0]     odd_bram1_b_rdata      ;
wire           [               7: 0]     even_bram2_b_rdata     ;
wire           [               7: 0]     odd_bram2_b_rdata      ;

bram_ture_dual_port
#(
        .C_ADDR_WIDTH                      (12                             ),
        .C_DATA_WIDTH                      (8                              ) 
)
u0_image_data_bram1
(
        .clka                              (clk_in1                        ),
        .wea                               (bram1_a_wenb                   ),
        .addra                             (bram_a_waddr                   ),
        .dina                              (bram_a_wdata                   ),
        .douta                             (                               ),
        .clkb                              (clk_in2                        ),
        .web                               (1'b0                           ),
        .addrb                             (even_bram1_b_raddr             ),
        .dinb                              (8'b0                           ),
        .doutb                             (even_bram1_b_rdata             ) 
);

bram_ture_dual_port
#(
        .C_ADDR_WIDTH                      (12                             ),
        .C_DATA_WIDTH                      (8                              ) 
)
u1_image_data_bram1
(
        .clka                              (clk_in1                        ),
        .wea                               (bram1_a_wenb                   ),
        .addra                             (bram_a_waddr                   ),
        .dina                              (bram_a_wdata                   ),
        .douta                             (                               ),
        .clkb                              (clk_in2                        ),
        .web                               (1'b0                           ),
        .addrb                             (odd_bram1_b_raddr              ),
        .dinb                              (8'b0                           ),
        .doutb                             (odd_bram1_b_rdata              ) 
);

bram_ture_dual_port
#(
        .C_ADDR_WIDTH                      (12                             ),
        .C_DATA_WIDTH                      (8                              ) 
)
u2_image_data_bram2
(
        .clka                              (clk_in1                        ),
        .wea                               (bram2_a_wenb                   ),
        .addra                             (bram_a_waddr                   ),
        .dina                              (bram_a_wdata                   ),
        .douta                             (                               ),
        .clkb                              (clk_in2                        ),
        .web                               (1'b0                           ),
        .addrb                             (even_bram2_b_raddr             ),
        .dinb                              (8'b0                           ),
        .doutb                             (even_bram2_b_rdata             ) 
);

bram_ture_dual_port
#(
        .C_ADDR_WIDTH                      (12                             ),
        .C_DATA_WIDTH                      (8                              ) 
)
u3_image_data_bram2
(
        .clka                              (clk_in1                        ),
        .wea                               (bram2_a_wenb                   ),
        .addra                             (bram_a_waddr                   ),
        .dina                              (bram_a_wdata                   ),
        .douta                             (                               ),
        .clkb                              (clk_in2                        ),
        .web                               (1'b0                           ),
        .addrb                             (odd_bram2_b_raddr              ),
        .dinb                              (8'b0                           ),
        .doutb                             (odd_bram2_b_rdata              ) 
);

wire                                     fifo_renb              ;
wire           [              10: 0]     fifo_rdata             ;
wire                                     fifo_empty             ;
wire                                     fifo_full              ;

asyn_fifo
#(
        .C_DATA_WIDTH                      (11                             ),
        .C_FIFO_DEPTH_WIDTH                (4                              ) 
)
u_tag_fifo
(
        .wr_rst_n                          (rst_n                          ),
        .wr_clk                            (clk_in1                        ),
        .wr_en                             (fifo_wenb                      ),
        .wr_data                           (fifo_wdata                     ),
        .wr_full                           (fifo_full                      ),
        .wr_cnt                            (                               ),
        .rd_rst_n                          (rst_n                          ),
        .rd_clk                            (clk_in2                        ),
        .rd_en                             (fifo_renb                      ),
        .rd_data                           (fifo_rdata                     ),
        .rd_empty                          (fifo_empty                     ),
        .rd_cnt                            (                               ) 
);

        localparam         S_IDLE           = 3'd0                 ;
        localparam         S_Y_LOAD         = 3'd1                 ;
        localparam         S_BRAM_ADDR      = 3'd2                 ;
        localparam         S_Y_INC          = 3'd3                 ;
        localparam         S_RD_FIFO        = 3'd4                 ;

reg            [               2: 0]     state                  ;
reg            [              26: 0]     y_dec                  ;
reg            [              26: 0]     x_dec                  ;
reg            [              10: 0]     y_cnt                  ;
reg            [              10: 0]     x_cnt                  ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
        state <= S_IDLE;
    else
    begin
        case(state)
            S_IDLE :
            begin
                if(fifo_empty == 1'b0)
                begin
                    if((fifo_rdata != 11'b0)&&(y_cnt == C_DST_IMG_HEIGHT))
                        state <= S_RD_FIFO;
                    else
                        state <= S_Y_LOAD;
                end
                else
                    state <= S_IDLE;
            end
            S_Y_LOAD :
            begin
                if((y_dec[26:16] + 1'b1 <= fifo_rdata)||(y_cnt == C_DST_IMG_HEIGHT - 1'b1))
                    state <= S_BRAM_ADDR;
                else
                    state <= S_RD_FIFO;
            end
            S_BRAM_ADDR :
            begin
                if(x_cnt == C_DST_IMG_WIDTH - 1'b1)
                    state <= S_Y_INC;
                else
                    state <= S_BRAM_ADDR;
            end
            S_Y_INC :
            begin
                if(y_cnt == C_DST_IMG_HEIGHT - 1'b1)
                    state <= S_RD_FIFO;
                else
                    state <= S_Y_LOAD;
            end
            S_RD_FIFO :
            begin
                state <= S_IDLE;
            end
            default :
            begin
                state <= S_IDLE;
            end
        endcase
    end
end

assign   fifo_renb        = (state == S_RD_FIFO) ? 1'b1 : 1'b0;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
        y_dec <= 27'b0;
    else
    begin
        if((state == S_IDLE)&&(fifo_empty == 1'b0)&&(fifo_rdata == 11'b0))
            y_dec <= 27'b0;
        else if(state == S_Y_INC)
            y_dec <= y_dec + C_Y_RATIO;
        else
            y_dec <= y_dec;
    end
end

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
        y_cnt <= 11'b0;
    else
    begin
        if((state == S_IDLE)&&(fifo_empty == 1'b0)&&(fifo_rdata == 11'b0))
            y_cnt <= 11'b0;
        else if(state == S_Y_INC)
            y_cnt <= y_cnt + 1'b1;
        else
            y_cnt <= y_cnt;
    end
end

always @(posedge clk_in2)
begin
    if(state == S_BRAM_ADDR)
        x_dec <= x_dec + C_X_RATIO;
    else
        x_dec <= 27'b0;
end

always @(posedge clk_in2)
begin
    if(state == S_BRAM_ADDR)
        x_cnt <= x_cnt + 1'b1;
    else
        x_cnt <= 11'b0;
end

//----------------------------------------------------------------------
//  c1
reg                                      img_vs_c1              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
        img_vs_c1 <= 1'b0;
    else
    begin
        if((state == S_BRAM_ADDR)&&(x_cnt == 11'b0)&&(y_cnt == 11'b0))
            img_vs_c1 <= 1'b1;
        else if((state == S_Y_INC)&&(y_cnt == C_DST_IMG_HEIGHT - 1'b1))
            img_vs_c1 <= 1'b0;
        else
            img_vs_c1 <= img_vs_c1;
    end
end

reg                                      img_hs_c1              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
        img_hs_c1 <= 1'b0;
    else
    begin
        if(state == S_BRAM_ADDR)
            img_hs_c1 <= 1'b1;
        else
            img_hs_c1 <= 1'b0;
    end
end

reg            [              10: 0]     x_int_c1               ;
reg            [              10: 0]     y_int_c1               ;
reg            [              16: 0]     x_fra_c1               ;
reg            [              16: 0]     inv_x_fra_c1           ;
reg            [              16: 0]     y_fra_c1               ;
reg            [              16: 0]     inv_y_fra_c1           ;

always @(posedge clk_in2)
begin
    x_int_c1     <= x_dec[25:16];
    y_int_c1     <= y_dec[25:16];
    x_fra_c1     <= {1'b0,x_dec[15:0]};
    inv_x_fra_c1 <= 17'h10000 - {1'b0,x_dec[15:0]};
    y_fra_c1     <= {1'b0,y_dec[15:0]};
    inv_y_fra_c1 <= 17'h10000 - {1'b0,y_dec[15:0]};
end

//----------------------------------------------------------------------
//  c2
reg                                      img_vs_c2              ;
reg                                      img_hs_c2              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c2 <= 1'b0;
        img_hs_c2 <= 1'b0;
    end
    else
    begin
        img_vs_c2 <= img_vs_c1;
        img_hs_c2 <= img_hs_c1;
    end
end

reg            [              11: 0]     bram_addr_c2           ;
reg            [              33: 0]     frac_00_c2             ;
reg            [              33: 0]     frac_01_c2             ;
reg            [              33: 0]     frac_10_c2             ;
reg            [              33: 0]     frac_11_c2             ;
reg                                      bram_mode_c2           ;

always @(posedge clk_in2)
begin
    bram_addr_c2 <= {y_int_c1[2:1],10'b0} + x_int_c1;
    frac_00_c2   <= inv_x_fra_c1 * inv_y_fra_c1;
    frac_01_c2   <= x_fra_c1 * inv_y_fra_c1;
    frac_10_c2   <= inv_x_fra_c1 * y_fra_c1;
    frac_11_c2   <= x_fra_c1 * y_fra_c1;
    bram_mode_c2 <= y_int_c1[0];
end

reg                                      right_pixel_extand_flag_c2  ;
reg                                      bottom_pixel_extand_flag_c2  ;

always @(posedge clk_in2)
begin
    if(x_int_c1 == C_SRC_IMG_WIDTH - 1'b1)
        right_pixel_extand_flag_c2 <= 1'b1;
    else
        right_pixel_extand_flag_c2 <= 1'b0;
    if(y_int_c1 == C_SRC_IMG_HEIGHT - 1'b1)
        bottom_pixel_extand_flag_c2 <= 1'b1;
    else
        bottom_pixel_extand_flag_c2 <= 1'b0;
end

//----------------------------------------------------------------------
//  c3
reg                                      img_vs_c3              ;
reg                                      img_hs_c3              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c3 <= 1'b0;
        img_hs_c3 <= 1'b0;
    end
    else
    begin
        img_vs_c3 <= img_vs_c2;
        img_hs_c3 <= img_hs_c2;
    end
end

always @(posedge clk_in2)
begin
    if(bram_mode_c2 == 1'b0)
    begin
        even_bram1_b_raddr <= bram_addr_c2;
        odd_bram1_b_raddr  <= bram_addr_c2 + 1'b1;
        even_bram2_b_raddr <= bram_addr_c2;
        odd_bram2_b_raddr  <= bram_addr_c2 + 1'b1;
    end
    else
    begin
        even_bram1_b_raddr <= bram_addr_c2 + 11'd1024;
        odd_bram1_b_raddr  <= bram_addr_c2 + 11'd1025;
        even_bram2_b_raddr <= bram_addr_c2;
        odd_bram2_b_raddr  <= bram_addr_c2 + 1'b1;
    end
end

reg            [              33: 0]     frac_00_c3             ;
reg            [              33: 0]     frac_01_c3             ;
reg            [              33: 0]     frac_10_c3             ;
reg            [              33: 0]     frac_11_c3             ;
reg                                      bram_mode_c3           ;
reg                                      right_pixel_extand_flag_c3  ;
reg                                      bottom_pixel_extand_flag_c3  ;

always @(posedge clk_in2)
begin
    frac_00_c3                  <= frac_00_c2;
    frac_01_c3                  <= frac_01_c2;
    frac_10_c3                  <= frac_10_c2;
    frac_11_c3                  <= frac_11_c2;
    bram_mode_c3                <= bram_mode_c2;
    right_pixel_extand_flag_c3  <= right_pixel_extand_flag_c2;
    bottom_pixel_extand_flag_c3 <= bottom_pixel_extand_flag_c2;
end

//----------------------------------------------------------------------
//  c4
reg                                      img_vs_c4              ;
reg                                      img_hs_c4              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c4 <= 1'b0;
        img_hs_c4 <= 1'b0;
    end
    else
    begin
        img_vs_c4 <= img_vs_c3;
        img_hs_c4 <= img_hs_c3;
    end
end

reg            [              33: 0]     frac_00_c4             ;
reg            [              33: 0]     frac_01_c4             ;
reg            [              33: 0]     frac_10_c4             ;
reg            [              33: 0]     frac_11_c4             ;
reg                                      bram_mode_c4           ;
reg                                      right_pixel_extand_flag_c4  ;
reg                                      bottom_pixel_extand_flag_c4  ;

always @(posedge clk_in2)
begin
    frac_00_c4                  <= frac_00_c3;
    frac_01_c4                  <= frac_01_c3;
    frac_10_c4                  <= frac_10_c3;
    frac_11_c4                  <= frac_11_c3;
    bram_mode_c4                <= bram_mode_c3;
    right_pixel_extand_flag_c4  <= right_pixel_extand_flag_c3;
    bottom_pixel_extand_flag_c4 <= bottom_pixel_extand_flag_c3;
end

//----------------------------------------------------------------------
//  c5
reg                                      img_vs_c5              ;
reg                                      img_hs_c5              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c5 <= 1'b0;
        img_hs_c5 <= 1'b0;
    end
    else
    begin
        img_vs_c5 <= img_vs_c4;
        img_hs_c5 <= img_hs_c4;
    end
end

reg            [               7: 0]     pixel_data00_c5        ;
reg            [               7: 0]     pixel_data01_c5        ;
reg            [               7: 0]     pixel_data10_c5        ;
reg            [               7: 0]     pixel_data11_c5        ;

always @(posedge clk_in2)
begin
    if(bram_mode_c4 == 1'b0)
    begin
        pixel_data00_c5 <= even_bram1_b_rdata;
        pixel_data01_c5 <= odd_bram1_b_rdata;
        pixel_data10_c5 <= even_bram2_b_rdata;
        pixel_data11_c5 <= odd_bram2_b_rdata;
    end
    else
    begin
        pixel_data00_c5 <= even_bram2_b_rdata;
        pixel_data01_c5 <= odd_bram2_b_rdata;
        pixel_data10_c5 <= even_bram1_b_rdata;
        pixel_data11_c5 <= odd_bram1_b_rdata;
    end
end

reg            [              33: 0]     frac_00_c5             ;
reg            [              33: 0]     frac_01_c5             ;
reg            [              33: 0]     frac_10_c5             ;
reg            [              33: 0]     frac_11_c5             ;
reg                                      right_pixel_extand_flag_c5  ;
reg                                      bottom_pixel_extand_flag_c5  ;

always @(posedge clk_in2)
begin
    frac_00_c5                  <= frac_00_c4;
    frac_01_c5                  <= frac_01_c4;
    frac_10_c5                  <= frac_10_c4;
    frac_11_c5                  <= frac_11_c4;
    right_pixel_extand_flag_c5  <= right_pixel_extand_flag_c4;
    bottom_pixel_extand_flag_c5 <= bottom_pixel_extand_flag_c4;
end

//----------------------------------------------------------------------
//  c6
reg                                      img_vs_c6              ;
reg                                      img_hs_c6              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c6 <= 1'b0;
        img_hs_c6 <= 1'b0;
    end
    else
    begin
        img_vs_c6 <= img_vs_c5;
        img_hs_c6 <= img_hs_c5;
    end
end

reg            [               7: 0]     pixel_data00_c6        ;
reg            [               7: 0]     pixel_data01_c6        ;
reg            [               7: 0]     pixel_data10_c6        ;
reg            [               7: 0]     pixel_data11_c6        ;

always @(posedge clk_in2)
begin
    case({right_pixel_extand_flag_c5,bottom_pixel_extand_flag_c5})
        2'b00 :
        begin
            pixel_data00_c6 <= pixel_data00_c5;
            pixel_data01_c6 <= pixel_data01_c5;
            pixel_data10_c6 <= pixel_data10_c5;
            pixel_data11_c6 <= pixel_data11_c5;
        end
        2'b01 :
        begin
            pixel_data00_c6 <= pixel_data00_c5;
            pixel_data01_c6 <= pixel_data01_c5;
            pixel_data10_c6 <= pixel_data00_c5;
            pixel_data11_c6 <= pixel_data01_c5;
        end
        2'b10 :
        begin
            pixel_data00_c6 <= pixel_data00_c5;
            pixel_data01_c6 <= pixel_data00_c5;
            pixel_data10_c6 <= pixel_data10_c5;
            pixel_data11_c6 <= pixel_data10_c5;
        end
        2'b11 :
        begin
            pixel_data00_c6 <= pixel_data00_c5;
            pixel_data01_c6 <= pixel_data00_c5;
            pixel_data10_c6 <= pixel_data00_c5;
            pixel_data11_c6 <= pixel_data00_c5;
        end
    endcase
end

reg            [              33: 0]     frac_00_c6             ;
reg            [              33: 0]     frac_01_c6             ;
reg            [              33: 0]     frac_10_c6             ;
reg            [              33: 0]     frac_11_c6             ;

always @(posedge clk_in2)
begin
    frac_00_c6 <= frac_00_c5;
    frac_01_c6 <= frac_01_c5;
    frac_10_c6 <= frac_10_c5;
    frac_11_c6 <= frac_11_c5;
end

//----------------------------------------------------------------------
//  c7
reg                                      img_vs_c7              ;
reg                                      img_hs_c7              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c7 <= 1'b0;
        img_hs_c7 <= 1'b0;
    end
    else
    begin
        img_vs_c7 <= img_vs_c6;
        img_hs_c7 <= img_hs_c6;
    end
end

reg            [              41: 0]     gray_data00_c7         ;
reg            [              41: 0]     gray_data01_c7         ;
reg            [              41: 0]     gray_data10_c7         ;
reg            [              41: 0]     gray_data11_c7         ;

always @(posedge clk_in2)
begin
    gray_data00_c7 <= frac_00_c6 * pixel_data00_c6;
    gray_data01_c7 <= frac_01_c6 * pixel_data01_c6;
    gray_data10_c7 <= frac_10_c6 * pixel_data10_c6;
    gray_data11_c7 <= frac_11_c6 * pixel_data11_c6;
end

//----------------------------------------------------------------------
//  c8
reg                                      img_vs_c8              ;
reg                                      img_hs_c8              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c8 <= 1'b0;
        img_hs_c8 <= 1'b0;
    end
    else
    begin
        img_vs_c8 <= img_vs_c7;
        img_hs_c8 <= img_hs_c7;
    end
end

reg            [              42: 0]     gray_data_tmp1_c8      ;
reg            [              42: 0]     gray_data_tmp2_c8      ;

always @(posedge clk_in2)
begin
    gray_data_tmp1_c8 <= gray_data00_c7 + gray_data01_c7;
    gray_data_tmp2_c8 <= gray_data10_c7 + gray_data11_c7;
end

//----------------------------------------------------------------------
//  c9
reg                                      img_vs_c9              ;
reg                                      img_hs_c9              ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c9 <= 1'b0;
        img_hs_c9 <= 1'b0;
    end
    else
    begin
        img_vs_c9 <= img_vs_c8;
        img_hs_c9 <= img_hs_c8;
    end
end

reg            [              43: 0]     gray_data_c9           ;

always @(posedge clk_in2)
begin
    gray_data_c9 <= gray_data_tmp1_c8 + gray_data_tmp2_c8;
end

//----------------------------------------------------------------------
//  c10
reg                                      img_vs_c10             ;
reg                                      img_hs_c10             ;

always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        img_vs_c10 <= 1'b0;
        img_hs_c10 <= 1'b0;
    end
    else
    begin
        img_vs_c10 <= img_vs_c9;
        img_hs_c10 <= img_hs_c9;
    end
end

reg            [              11: 0]     gray_data_c10          ;

always @(posedge clk_in2)
begin
    gray_data_c10 <= gray_data_c9[43:32] + gray_data_c9[31];
end

//----------------------------------------------------------------------
//  signals output
always @(posedge clk_in2)
begin
    if(rst_n == 1'b0)
    begin
        post_img_vsync <= 1'b0;
        post_img_href  <= 1'b0;
    end
    else
    begin
        post_img_vsync <= img_vs_c10;
        post_img_href  <= img_hs_c10;
    end
end

always @(posedge clk_in2)
begin
    if(gray_data_c10 > 12'd255)
        post_img_gray <= 8'd255;
    else
        post_img_gray <= gray_data_c10[7:0];
end

endmodule
```

#### 缓存buff设计
![[ISP之常用图像缩放算法 2025-08-03 19.57.35.excalidraw|1200]]
#### 效果
![image.png](https://obsidian-picturebed-1256135654.cos.ap-nanjing.myqcloud.com/obsidion/20250807111709516.png?imageSlim)


# BUG日志
*“今天又在很努力的写bug了！”*
#### 1. 在Verilog实现地址生成的过程中，难以避免对于除法精度的处理 
如下公式，往往会出现小数，就算对其扩大，由于dst也是较大的数字，一旦误差被放大1000也是巨大的偏移，会越来越偏，花了大半天时间，很难在资源和效率之间平衡，每次都用除法资源量是巨大的，不可接受！！又想着去修正数值，但是这样又会写死，难以适配更多的变化。

最后，取消了`-1` 🤤，不可否认的是这样也存在难以适配的问题，但是这也许是当先最简单最省事最可靠的办法，对于其它分辨率变化，我认为更应该是在 `localparam` 部分就让他为整数。

$$
\begin{align*}
src_x &= dst_x \cdot \frac{W_s - 1}{W_d - 1} \\
src_y &= dst_y \cdot \frac{H_s - 1}{H_d - 1}
\end{align*}
$$
#### 2. 反压机制存在问题，在消隐期间启动
应该需要一个等待机制，用来匹配速度。

**机制原理：** 当`addr_gen`需要的数据还没有被输入视频流准备好的时候，缓冲（`line_buffer`）必须通知`addr_gen`，`addr_gen`则必须暂停，直到数据就绪、

