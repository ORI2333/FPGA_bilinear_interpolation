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
        parameter                          MODE             = "GRAYSCALE",  // ����ģʽ: "RGB" �� "GRAYSCALE"
        parameter                          IMG_WIDTH        = 640   ,        
        parameter                          IMG_HEIGHT       = 480   ,        
        parameter                          PIXEL_WIDTH      = 8     ,         
        parameter                          FRAMES_TO_DUMP   = 1     ,          // Ҫת������֡��
        parameter                          OUTPUT_FILE_BASE = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/FPGA_output_hex"// ��������ļ���
) (

        input  wire                                       clk                 ,// ϵͳʱ��
        input  wire                                       rst_n               ,// �͵�ƽ��Ч�첽��λ

        input  wire                                       i_vsync             ,// ��ͬ���ź� (��ǰδʹ�ã�������)
        input  wire                                       i_hsync             ,// ��ͬ���ź� (��ǰδʹ�ã�������)
        input  wire                                       i_data_valid        ,// ��Ч�������ݱ�־
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_r            ,// R���� �� �Ҷ�ֵ
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_g            ,// G����
        input  wire         [   PIXEL_WIDTH-1: 0]         i_data_b            ,// B����

        output reg          [              31: 0]         o_frame_count        

);


        localparam                         TOTAL_PIXELS     = IMG_WIDTH * IMG_HEIGHT;

// --- ״̬�� ---
        localparam                         S_IDLE           = 1'b0  ; // ����״̬
        localparam                         S_DUMPING        = 1'b1  ; // ����ת��״̬
reg                                      state                  ;

// --- �ļ����� ---
integer                                  file_handle            ;
reg            [         8*128-1: 0]     current_filename       ;// ���ڴ洢��̬���ɵ��ļ���

// --- ������ ---
reg            [              31: 0]     pixel_count            ;// ��д������ؼ�����
reg            [              31: 0]     frame_count            ;// ����ɵ�֡������

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
                // �ڿ���״̬���ȴ���һ����Ч����
                if (i_data_valid && (frame_count < FRAMES_TO_DUMP)) begin
                    // --- ״̬�л� & ���ļ� ---
                    state <= S_DUMPING;

                    // ��̬�����ļ���
                    $sformat(current_filename, "%s_%0d.txt", OUTPUT_FILE_BASE, frame_count);
                    file_handle = $fopen(current_filename, "w");
                    if (file_handle == 0) begin
                        $display("Error: Dumper failed to open %s", current_filename);
                        $stop;
                    end else begin
                        $display("Dumper: Frame %0d. Output file '%s' opened.", frame_count, current_filename);
                    end

                    // --- ����д���һ������ ---
                    if (MODE == "RGB") begin
                        $fdisplay(file_handle, "%h %h %h", i_data_r, i_data_g, i_data_b);
                    end else begin                                  // GRAYSCALE
                        $fdisplay(file_handle, "%h", i_data_r);
                    end
                    pixel_count <= 1;                               // ���ؼ�������1��ʼ
                end
            end

            S_DUMPING: begin
                // ��ת��״̬������д����Ч����
                if (i_data_valid) begin
                    if (MODE == "RGB") begin
                        $fdisplay(file_handle, "%h %h %h", i_data_r, i_data_g, i_data_b);
                    end else begin                                  // GRAYSCALE
                        $fdisplay(file_handle, "%h", i_data_r);
                    end

                    // ����Ƿ���д��һ֡
                    if (pixel_count == TOTAL_PIXELS - 1) begin
                        // --- ״̬�л� & �ر��ļ� ---
                        state <= S_IDLE;
                        $display("Dumper: Frame %0d dumping complete. Closing file.", frame_count);
                        $fclose(file_handle);
                        file_handle = 0;
                        frame_count <= frame_count + 1;             // ���һ֡��֡��������1
                        pixel_count <= 0;
                    end else begin
                        pixel_count <= pixel_count + 1;
                    end
                end
            end
        endcase
    end
end

assign o_frame_count = frame_count; // �����ǰ֡����

endmodule

/*
// tb.sv
// ������
// ��������
always @(posedge clk) begin
    if (o_frame_count >= FRAMES_TO_DUMP) begin
        $display("All frames have been dumped (%0d frames), simulation finished.", FRAMES_TO_DUMP);
        $finish;
    end
end
*/
