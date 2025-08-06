`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 2025/07/26 20:17:52
// Design Name: 
// Module Name: bilinear_scaler_top
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 从TXT文件读取像素数据，生成模拟相机视频流，用于仿真测试。
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////
module image_gen #(
        parameter                          MODE             = "GRAYSCALE",           // 工作模式: "RGB" 或 "GRAYSCALE"
        parameter                          IMG_WIDTH        = 640   ,           // 图像有效区域宽度
        parameter                          IMG_HEIGHT       = 480   ,           // 图像有效区域高度
        parameter                          H_BLANK          = 160   ,           // 行消隐周期 (包括 hsync 脉冲宽度)
        parameter                          V_BLANK          = 45    ,           // 场消隐周期 (包括 vsync 脉冲宽度)
 //       parameter                          H_SYNC_PULSE     = 96    ,           // hsync 脉冲宽度
 //       parameter                          V_SYNC_PULSE     = 2     ,           // vsync 脉冲宽度
        parameter                          PIXEL_WIDTH      = 8     ,           // 每个颜色分量的位宽
        parameter                          DATA_FILE        = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/output_pixels_hex.txt"// 输入的像素数据文件名
)(

        input  wire                                       clk                 ,// 系统时钟
        input  wire                                       rst_n               ,// 低电平有效异步复位

        output reg                                        o_vsync             ,// 场同步信号
        output reg                                        o_hsync             ,// 行同步信号
        output reg                                        o_data_valid        ,// 有效像素数据标志
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_r            ,// R分量 或 灰度值
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_g            ,// G分量 (灰度模式下与R相同)
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_b             // B分量 (灰度模式下与R相同)
);

// =================================================================
// 内部信号与常量定义 (Internal Signals & Constants)
// =================================================================
        localparam                         H_TOTAL          = IMG_WIDTH + H_BLANK; // 一行的总周期数
        localparam                         V_TOTAL          = IMG_HEIGHT + V_BLANK; // 一帧的总行数
        localparam                         TOTAL_PIXELS     = IMG_WIDTH * IMG_HEIGHT;

// 行、场计数器
reg            [              15: 0]     h_count                ;
reg            [              15: 0]     v_count                ;

// 像素数据存储器
reg            [   PIXEL_WIDTH-1: 0]     mem_r[0:TOTAL_PIXELS-1]  ;
reg            [   PIXEL_WIDTH-1: 0]     mem_g[0:TOTAL_PIXELS-1]  ;
reg            [   PIXEL_WIDTH-1: 0]     mem_b[0:TOTAL_PIXELS-1]  ;

// =================================================================
// 文件读取逻辑 (File Reading Logic)
// =================================================================
initial begin
    $display("Stimulus Generator: Loading pixel data from '%s' in %s mode.", DATA_FILE, MODE);
    if (MODE == "RGB") begin
        // 读取RGB三通道数据
        $readmemh(DATA_FILE, mem_r, 0, TOTAL_PIXELS-1);             // 使用 $readmemh 读取更高效
        // 注意: $readmemh 无法直接读取多列, 这里假设文件格式为 R G B 交错
        // 更稳健的方式是使用 $fscanf 循环读取
        // integer file_id, i;
        // file_id = $fopen(DATA_FILE, "r");
        // if (file_id == 0) begin
        //     $display("Error: Failed to open %s", DATA_FILE); $stop;
        // end
        // for (i = 0; i < TOTAL_PIXELS; i = i + 1) begin
        //     $fscanf(file_id, "%d %d %d\n", mem_r[i], mem_g[i], mem_b[i]);
        // end
        // $fclose(file_id);
    end else if (MODE == "GRAYSCALE") begin
        // 读取单通道灰度数据
        $readmemh(DATA_FILE, mem_r);                                // 对于灰度图，只填充R通道内存
    end
    $display("Stimulus Generator: Pixel data loaded successfully.");
end


// 时序生成逻辑
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        h_count <= 0;
        v_count <= 0;
    end else begin
        // 行计数器
        if (h_count < H_TOTAL - 1) begin
            h_count <= h_count + 1;
        end else begin
            h_count <= 0;
            // 场计数器
            if (v_count < V_TOTAL - 1) begin
                v_count <= v_count + 1;
            end else
            begin
                v_count <= 0;
            end
        end
    end
end




reg            [              15: 0]     active_h_count         ;
reg            [              15: 0]     active_v_count         ;
wire                                     next_data_valid        ;
// 有效区域计数器，用于正确寻址
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        active_h_count <= 0;
        active_v_count <= 0;
    end else begin
        // 在一帧开始前，复位active计数器
        if (v_count == V_TOTAL - 1 && h_count == H_TOTAL - 1) begin
             active_h_count <= 0;
             active_v_count <= 0;
        end else if (next_data_valid) begin                         // 只在有效数据期间计数
            if (active_h_count < IMG_WIDTH - 1) begin
                active_h_count <= active_h_count + 1;
            end else begin
                active_h_count <= 0;
                if (active_v_count < IMG_HEIGHT - 1) begin
                    active_v_count <= active_v_count + 1;
                end else begin
                    active_v_count <= 0;
                end
            end
        end
    end
end
assign         next_data_valid  = (h_count < IMG_WIDTH) && (v_count < IMG_HEIGHT);


wire           [              31: 0]     pixel_index          =active_v_count * IMG_WIDTH + active_h_count;

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        o_vsync      <= 1'b1;                                       // 复位到消隐状态
        o_hsync      <= 1'b1;
        o_data_valid <= 1'b0;
        o_data_r     <= 'd0;
        o_data_g     <= 'd0;
        o_data_b     <= 'd0;
    end else begin
        // 寄存所有输出，确保与时钟同步
        o_hsync <= (h_count >= IMG_WIDTH);
        o_vsync <= (v_count >= IMG_HEIGHT);
        o_data_valid <= next_data_valid;

        if (next_data_valid) begin
            if (MODE == "RGB") begin
                o_data_r <= mem_r[pixel_index];
                o_data_g <= mem_g[pixel_index];
                o_data_b <= mem_b[pixel_index];
            end else begin                                          // GRAYSCALE mode
                o_data_r <= mem_r[pixel_index];
                o_data_g <= mem_r[pixel_index];
                o_data_b <= mem_r[pixel_index];
            end
        end else begin
             o_data_r <= 'd0;
             o_data_g <= 'd0;
             o_data_b <= 'd0;
        end
    end
end


endmodule