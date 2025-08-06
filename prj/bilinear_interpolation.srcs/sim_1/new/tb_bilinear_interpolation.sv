`timescale 1ns / 1ps

module tb_bilinear_interpolation;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 参数定义 (Parameters)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 源图像参数
        localparam         SRC_IMG_WIDTH    = 640                  ;  // 源图像宽度
        localparam         SRC_IMG_HEIGHT   = 480                  ;  // 源图像高度
        localparam         SRC_H_BLANK      = 160                  ;  // 源图像行消隐
        localparam         SRC_V_BLANK      = 45                   ;   // 源图像场消隐

    // 目标图像参数
        localparam         DST_IMG_WIDTH    = 1024                 ; // 目标图像宽度
        localparam         DST_IMG_HEIGHT   = 768                  ;  // 目标图像高度

    // 缩放比例 (根据DUT注释预先计算)
    // floor(SRC_W / DST_W * 2^16) = floor(640/1024 * 65536) = 40960
        localparam         X_RATIO          = 16'd40960            ;
    // floor(SRC_H / DST_H * 2^16) = floor(480/768 * 65536) = 40960
        localparam         Y_RATIO          = 16'd40960            ;

    // 数据位宽
        localparam         PIXEL_WIDTH      = 8                    ;

    // 仿真控制
        localparam         FRAMES_TO_SIM    = 1                    ;    // 设置要仿真的总帧数
    localparam string INPUT_FILE_NAME  = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/output_pixels_hex.txt";   // 输入文件名
    localparam string OUTPUT_FILE_BASE = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/FPGA_output_hex.txt"; // 输出文件名前缀

    // 时钟周期
        localparam         CLK1_PERIOD      = 20.0                 ; // 50 MHz, 用于输入级 (clk_in1)
        localparam         CLK2_PERIOD      = 10.0                 ; // 100 MHz, 用于处理/输出级 (clk_in2)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 信号定义 (Signals)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reg                                      clk_in1                ;// 时钟1
reg                                      clk_in2                ;// 时钟2
reg                                      rst_n                  ;// 复位信号 (低有效)

    // 从 image_gen 到 DUT 的信号
wire                                     src_vsync_raw          ;// 来自 image_gen, 消隐期为高
wire                                     src_hsync_raw          ;// 来自 image_gen, 消隐期为高
wire                                     src_href               ;// 来自 image_gen 的数据有效信号
wire           [   PIXEL_WIDTH-1: 0]     src_gray               ;// 源灰度数据

    // 从 DUT 到 image_return 的信号
wire                                     dst_vsync              ;// 目标场同步
wire                                     dst_href               ;// 目标行有效
wire           [   PIXEL_WIDTH-1: 0]     dst_gray               ;// 目标灰度数据

    // 从 image_return 返回的仿真状态信号
wire           [              31: 0]     frames_dumped_count    ;// 已转储的帧数


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 时钟和复位生成 (Clock and Reset Generation)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    initial begin
        clk_in1 = 1'b0;
        forever #(CLK1_PERIOD/2) clk_in1 = ~clk_in1;
    end

    initial begin
        clk_in2 = 1'b0;
        forever #(CLK2_PERIOD/2) clk_in2 = ~clk_in2;
    end

    initial begin
        $display("-------------------------------------------");
        $display("--- 开始双线性插值 Testbench ---");
        $display("-------------------------------------------");

        // 自动创建测试图像文件
       // generate_input_image(INPUT_FILE_NAME, SRC_IMG_WIDTH, SRC_IMG_HEIGHT);

        // 施加复位
        rst_n = 1'b0;
        #(CLK1_PERIOD * 10);
        rst_n = 1'b1;
        $display("[%0t ns] 复位已释放。", $time);
    end

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 模块实例化 (Module Instantiation)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 激励生成器 (从文件读取数据并创建视频流)
    image_gen #(
        .MODE                              ("GRAYSCALE"                    ),
        .IMG_WIDTH                         (SRC_IMG_WIDTH                  ),
        .IMG_HEIGHT                        (SRC_IMG_HEIGHT                 ),
        .H_BLANK                           (SRC_H_BLANK                    ),
        .V_BLANK                           (SRC_V_BLANK                    ),
        .PIXEL_WIDTH                       (PIXEL_WIDTH                    ),
        .DATA_FILE                         (INPUT_FILE_NAME                ) 
    ) u_image_gen (
        .clk                               (clk_in1                        ),
        .rst_n                             (rst_n                          ),
        .o_vsync                           (src_vsync_raw                  ),// 注意: 此信号在消隐期间为高
        .o_hsync                           (src_hsync_raw                  ),
        .o_data_valid                      (src_href                       ),
        .o_data_r                          (src_gray                       ),
        .o_data_g                          (                               ),// 灰度模式下未使用
        .o_data_b                          (                               ) // 灰度模式下未使用
    );
    
    // 被测设计 (Design Under Test, DUT)
    bilinear_interpolation #(
        .C_SRC_IMG_WIDTH                   (SRC_IMG_WIDTH                  ),
        .C_SRC_IMG_HEIGHT                  (SRC_IMG_HEIGHT                 ),
        .C_DST_IMG_WIDTH                   (DST_IMG_WIDTH                  ),
        .C_DST_IMG_HEIGHT                  (DST_IMG_HEIGHT                 ),
        .C_X_RATIO                         (X_RATIO                        ),
        .C_Y_RATIO                         (Y_RATIO                        ) 
    ) u_dut (
        .clk_in1                           (clk_in1                        ),
        .clk_in2                           (clk_in2                        ),
        .rst_n                             (rst_n                          ),

        // 输入图像流 (连接到激励生成器)
        // DUT的 per_img_vsync 需要在有效帧期间为高, 所以我们反转 image_gen 的输出
        .per_img_vsync                     (~src_vsync_raw                 ),
        .per_img_href                      (src_href                       ),
        .per_img_gray                      (src_gray                       ),
        
        // 输出图像流 (连接到响应检查器)
        .post_img_vsync                    (dst_vsync                      ),
        .post_img_href                     (dst_href                       ),
        .post_img_gray                     (dst_gray                       ) 
    );

    // 响应检查器/转储器 (接收视频流并写入文件)
    image_return #(
        .MODE                              ("GRAYSCALE"                    ),
        .IMG_WIDTH                         (DST_IMG_WIDTH                  ),
        .IMG_HEIGHT                        (DST_IMG_HEIGHT                 ),
        .PIXEL_WIDTH                       (PIXEL_WIDTH                    ),
        .FRAMES_TO_DUMP                    (FRAMES_TO_SIM                  ),
        .OUTPUT_FILE_BASE                  (OUTPUT_FILE_BASE               ) 
    ) u_image_return (
        .clk                               (clk_in2                        ),
        .rst_n                             (rst_n                          ),
        .i_vsync                           (dst_vsync                      ),
        .i_hsync                           (~dst_href                      ),// 连接到一个合理的值, 不关键
        .i_data_valid                      (dst_href                       ),
        .i_data_r                          (dst_gray                       ),
        .i_data_g                          (dst_gray                       ),
        .i_data_b                          (dst_gray                       ),
        .o_frame_count                     (frames_dumped_count            ) 
    );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // 仿真控制与任务 (Simulation Control & Tasks)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 在所需数量的帧被转储后停止仿真
    always @(posedge clk_in2) begin
        if (frames_dumped_count >= FRAMES_TO_SIM) begin
            $display("-------------------------------------------");
            $display("[%0t ns] 仿真成功结束。", $time);
            $display("--- 已转储 %0d 帧图像。 ---", FRAMES_TO_SIM);
            $display("--- 输出文件: %s_%0d.txt ---", OUTPUT_FILE_BASE, FRAMES_TO_SIM-1);
            $display("-------------------------------------------");
            $finish;
        end
    end
    
    // 添加一个超时机制, 防止仿真因意外错误而永远运行下去
    initial begin
        #(CLK1_PERIOD * (SRC_IMG_WIDTH + SRC_H_BLANK) * (SRC_IMG_HEIGHT + SRC_V_BLANK) * (FRAMES_TO_SIM + 1));
        $display("-------------------------------------------");
        $error("[%0t ns] 仿真超时！仿真运行时间过长。", $time);
        $display("-------------------------------------------");
        $finish;
    end

    // 用于生成测试图像文件 (例如, 一个水平渐变图像) 的任务
    task automatic generate_input_image(string filename, integer width, integer height);
integer                                  file_handle            ;
integer                                  i,                   j;
integer                                  gray_value             ;
        
        file_handle = $fopen(filename, "w");
        if (file_handle == 0) begin
            $error("TB 错误: 无法打开文件 '%s' 进行写入。", filename);
            $finish;
        end

        $display("TB 信息: 正在生成测试图像 '%s' (%0d x %0d)...", filename, width, height);
        for (i = 0; i < height; i = i + 1) begin
            for (j = 0; j < width; j = j + 1) begin
                // 一个从0到255的简单水平渐变
                gray_value = (j * 255) / (width - 1);
                $fwrite(file_handle, "%x\n", gray_value);
            end
        end
        $fclose(file_handle);
        $display("TB 信息: 测试图像生成完毕。");
    endtask

endmodule