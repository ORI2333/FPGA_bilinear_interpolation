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
// ��TXT�ļ���ȡ�������ݣ�����ģ�������Ƶ�������ڷ�����ԡ�
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////
module image_gen #(
        parameter                          MODE             = "GRAYSCALE",           // ����ģʽ: "RGB" �� "GRAYSCALE"
        parameter                          IMG_WIDTH        = 640   ,           // ͼ����Ч������
        parameter                          IMG_HEIGHT       = 480   ,           // ͼ����Ч����߶�
        parameter                          H_BLANK          = 160   ,           // ���������� (���� hsync ������)
        parameter                          V_BLANK          = 45    ,           // ���������� (���� vsync ������)
 //       parameter                          H_SYNC_PULSE     = 96    ,           // hsync ������
 //       parameter                          V_SYNC_PULSE     = 2     ,           // vsync ������
        parameter                          PIXEL_WIDTH      = 8     ,           // ÿ����ɫ������λ��
        parameter                          DATA_FILE        = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/output_pixels_hex.txt"// ��������������ļ���
)(

        input  wire                                       clk                 ,// ϵͳʱ��
        input  wire                                       rst_n               ,// �͵�ƽ��Ч�첽��λ

        output reg                                        o_vsync             ,// ��ͬ���ź�
        output reg                                        o_hsync             ,// ��ͬ���ź�
        output reg                                        o_data_valid        ,// ��Ч�������ݱ�־
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_r            ,// R���� �� �Ҷ�ֵ
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_g            ,// G���� (�Ҷ�ģʽ����R��ͬ)
        output reg          [   PIXEL_WIDTH-1: 0]         o_data_b             // B���� (�Ҷ�ģʽ����R��ͬ)
);

// =================================================================
// �ڲ��ź��볣������ (Internal Signals & Constants)
// =================================================================
        localparam                         H_TOTAL          = IMG_WIDTH + H_BLANK; // һ�е���������
        localparam                         V_TOTAL          = IMG_HEIGHT + V_BLANK; // һ֡��������
        localparam                         TOTAL_PIXELS     = IMG_WIDTH * IMG_HEIGHT;

// �С���������
reg            [              15: 0]     h_count                ;
reg            [              15: 0]     v_count                ;

// �������ݴ洢��
reg            [   PIXEL_WIDTH-1: 0]     mem_r[0:TOTAL_PIXELS-1]  ;
reg            [   PIXEL_WIDTH-1: 0]     mem_g[0:TOTAL_PIXELS-1]  ;
reg            [   PIXEL_WIDTH-1: 0]     mem_b[0:TOTAL_PIXELS-1]  ;

// =================================================================
// �ļ���ȡ�߼� (File Reading Logic)
// =================================================================
initial begin
    $display("Stimulus Generator: Loading pixel data from '%s' in %s mode.", DATA_FILE, MODE);
    if (MODE == "RGB") begin
        // ��ȡRGB��ͨ������
        $readmemh(DATA_FILE, mem_r, 0, TOTAL_PIXELS-1);             // ʹ�� $readmemh ��ȡ����Ч
        // ע��: $readmemh �޷�ֱ�Ӷ�ȡ����, ��������ļ���ʽΪ R G B ����
        // ���Ƚ��ķ�ʽ��ʹ�� $fscanf ѭ����ȡ
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
        // ��ȡ��ͨ���Ҷ�����
        $readmemh(DATA_FILE, mem_r);                                // ���ڻҶ�ͼ��ֻ���Rͨ���ڴ�
    end
    $display("Stimulus Generator: Pixel data loaded successfully.");
end


// ʱ�������߼�
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        h_count <= 0;
        v_count <= 0;
    end else begin
        // �м�����
        if (h_count < H_TOTAL - 1) begin
            h_count <= h_count + 1;
        end else begin
            h_count <= 0;
            // ��������
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
// ��Ч�����������������ȷѰַ
always @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
        active_h_count <= 0;
        active_v_count <= 0;
    end else begin
        // ��һ֡��ʼǰ����λactive������
        if (v_count == V_TOTAL - 1 && h_count == H_TOTAL - 1) begin
             active_h_count <= 0;
             active_v_count <= 0;
        end else if (next_data_valid) begin                         // ֻ����Ч�����ڼ����
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
        o_vsync      <= 1'b1;                                       // ��λ������״̬
        o_hsync      <= 1'b1;
        o_data_valid <= 1'b0;
        o_data_r     <= 'd0;
        o_data_g     <= 'd0;
        o_data_b     <= 'd0;
    end else begin
        // �Ĵ����������ȷ����ʱ��ͬ��
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