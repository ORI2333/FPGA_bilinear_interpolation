`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 2025/07/27 18:52:47
// Design Name: 
// Module Name: image_return
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module image_return #(
        parameter                          MODE             = "GRAYSCALE",  // 工作模式: "RGB" 或 "GRAYSCALE"
        parameter                          IMG_WIDTH        = 640   ,        
        parameter                          IMG_HEIGHT       = 480   ,        
        parameter                          PIXEL_WIDTH      = 8     ,         
        parameter                          FRAMES_TO_DUMP   = 1     ,          // 要转储的总帧数
        parameter                          OUTPUT_FILE_BASE = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/FPGA_output_hex"// 输出数据文件名
) (

        input  wire                                       clk                 ,// 系统时钟
        input  wire                                       rst_n               ,// 低电平有效异步复位

        input  wire                                       i_vsync             ,// 场同步信号 (当前未使用，但保留)
        input  wire                                       i_hsync             ,// 行同步信号 (当前未使用，但保留)
        input  wire                                       i_data_valid        ,// 有效像素数据标志
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_r            ,// R分量 或 灰度值
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_g            ,// G分量
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_b            ,// B分量

        output reg          [              31: 0]         o_frame_count        

);


        localparam                         TOTAL_PIXELS     = IMG_WIDTH * IMG_HEIGHT;

// --- 状态机 ---
        localparam                         S_IDLE           = 1'b0  ; // 空闲状态
        localparam                         S_DUMPING        = 1'b1  ; // 数据转储状态
reg                                      state                  ;

// --- 文件操作 ---
integer                                  file_handle            ;
reg            [         8*128-1: 0]     current_filename       ;// 用于存储动态生成的文件名

// --- 计数器 ---
reg            [              31: 0]     pixel_count            ;// 已写入的像素计数器
reg            [              31: 0]     frame_count            ;// 已完成的帧计数器

always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        state       <= S_IDLE;
        pixel_count <= 0;
        frame_count <= 0;
        if (file_handle) begin
            $fclose(file_handle);
            file_handle = 0;
        end
    end else begin
        case (state)
            S_IDLE: begin
                // 在空闲状态，等待第一个有效数据
                if (i_data_valid && (frame_count < FRAMES_TO_DUMP)) begin
                    // --- 状态切换 & 打开文件 ---
                    state <= S_DUMPING;

                    // 动态生成文件名
                    $sformat(current_filename, "%s_%0d.txt", OUTPUT_FILE_BASE, frame_count);
                    file_handle = $fopen(current_filename, "w");
                    if (file_handle == 0) begin
                        $display("Error: Dumper failed to open %s", current_filename);
                        $stop;
                    end else begin
                        $display("Dumper: Frame %0d. Output file '%s' opened.", frame_count, current_filename);
                    end

                    // --- 立即写入第一个像素 ---
                    if (MODE == "RGB") begin
                        $fdisplay(file_handle, "%h %h %h", i_data_r, i_data_g, i_data_b);
                    end else begin                                  // GRAYSCALE
                        $fdisplay(file_handle, "%h", i_data_r);
                    end
                    pixel_count <= 1;                               // 像素计数器从1开始
                end
            end

            S_DUMPING: begin
                // 在转储状态，继续写入有效数据
                if (i_data_valid) begin
                    if (MODE == "RGB") begin
                        $fdisplay(file_handle, "%h %h %h", i_data_r, i_data_g, i_data_b);
                    end else begin                                  // GRAYSCALE
                        $fdisplay(file_handle, "%h", i_data_r);
                    end

                    // 检查是否已写完一帧
                    if (pixel_count == TOTAL_PIXELS - 1) begin
                        // --- 状态切换 & 关闭文件 ---
                        state <= S_IDLE;
                        $display("Dumper: Frame %0d dumping complete. Closing file.", frame_count);
                        $fclose(file_handle);
                        file_handle = 0;
                        frame_count <= frame_count + 1;             // 完成一帧，帧计数器加1
                        pixel_count <= 0;
                    end else begin
                        pixel_count <= pixel_count + 1;
                    end
                end
            end
        endcase
    end
end

assign o_frame_count = frame_count; // 输出当前帧计数

endmodule

/*
// tb.sv
// 测试用
// 结束条件
always @(posedge clk) begin
    if (o_frame_count >= FRAMES_TO_DUMP) begin
        $display("All frames have been dumped (%0d frames), simulation finished.", FRAMES_TO_DUMP);
        $finish;
    end
end
*/
