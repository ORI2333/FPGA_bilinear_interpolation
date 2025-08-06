`timescale 1ns / 1ps

module tb_bilinear_interpolation;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // �������� (Parameters)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Դͼ�����
        localparam         SRC_IMG_WIDTH    = 640                  ;  // Դͼ����
        localparam         SRC_IMG_HEIGHT   = 480                  ;  // Դͼ��߶�
        localparam         SRC_H_BLANK      = 160                  ;  // Դͼ��������
        localparam         SRC_V_BLANK      = 45                   ;   // Դͼ������

    // Ŀ��ͼ�����
        localparam         DST_IMG_WIDTH    = 1024                 ; // Ŀ��ͼ����
        localparam         DST_IMG_HEIGHT   = 768                  ;  // Ŀ��ͼ��߶�

    // ���ű��� (����DUTע��Ԥ�ȼ���)
    // floor(SRC_W / DST_W * 2^16) = floor(640/1024 * 65536) = 40960
        localparam         X_RATIO          = 16'd40960            ;
    // floor(SRC_H / DST_H * 2^16) = floor(480/768 * 65536) = 40960
        localparam         Y_RATIO          = 16'd40960            ;

    // ����λ��
        localparam         PIXEL_WIDTH      = 8                    ;

    // �������
        localparam         FRAMES_TO_SIM    = 1                    ;    // ����Ҫ�������֡��
    localparam string INPUT_FILE_NAME  = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/output_pixels_hex.txt";   // �����ļ���
    localparam string OUTPUT_FILE_BASE = "F:/EngineeringWarehouse/ISP/ImageGen/ImageGenPY/FPGA_output_hex.txt"; // ����ļ���ǰ׺

    // ʱ������
        localparam         CLK1_PERIOD      = 20.0                 ; // 50 MHz, �������뼶 (clk_in1)
        localparam         CLK2_PERIOD      = 10.0                 ; // 100 MHz, ���ڴ���/����� (clk_in2)

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // �źŶ��� (Signals)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
reg                                      clk_in1                ;// ʱ��1
reg                                      clk_in2                ;// ʱ��2
reg                                      rst_n                  ;// ��λ�ź� (����Ч)

    // �� image_gen �� DUT ���ź�
wire                                     src_vsync_raw          ;// ���� image_gen, ������Ϊ��
wire                                     src_hsync_raw          ;// ���� image_gen, ������Ϊ��
wire                                     src_href               ;// ���� image_gen ��������Ч�ź�
wire           [   PIXEL_WIDTH-1: 0]     src_gray               ;// Դ�Ҷ�����

    // �� DUT �� image_return ���ź�
wire                                     dst_vsync              ;// Ŀ�곡ͬ��
wire                                     dst_href               ;// Ŀ������Ч
wire           [   PIXEL_WIDTH-1: 0]     dst_gray               ;// Ŀ��Ҷ�����

    // �� image_return ���صķ���״̬�ź�
wire           [              31: 0]     frames_dumped_count    ;// ��ת����֡��


    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ʱ�Ӻ͸�λ���� (Clock and Reset Generation)
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
        $display("--- ��ʼ˫���Բ�ֵ Testbench ---");
        $display("-------------------------------------------");

        // �Զ���������ͼ���ļ�
       // generate_input_image(INPUT_FILE_NAME, SRC_IMG_WIDTH, SRC_IMG_HEIGHT);

        // ʩ�Ӹ�λ
        rst_n = 1'b0;
        #(CLK1_PERIOD * 10);
        rst_n = 1'b1;
        $display("[%0t ns] ��λ���ͷš�", $time);
    end

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ģ��ʵ���� (Module Instantiation)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ���������� (���ļ���ȡ���ݲ�������Ƶ��)
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
        .o_vsync                           (src_vsync_raw                  ),// ע��: ���ź��������ڼ�Ϊ��
        .o_hsync                           (src_hsync_raw                  ),
        .o_data_valid                      (src_href                       ),
        .o_data_r                          (src_gray                       ),
        .o_data_g                          (                               ),// �Ҷ�ģʽ��δʹ��
        .o_data_b                          (                               ) // �Ҷ�ģʽ��δʹ��
    );
    
    // ������� (Design Under Test, DUT)
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

        // ����ͼ���� (���ӵ�����������)
        // DUT�� per_img_vsync ��Ҫ����Ч֡�ڼ�Ϊ��, �������Ƿ�ת image_gen �����
        .per_img_vsync                     (~src_vsync_raw                 ),
        .per_img_href                      (src_href                       ),
        .per_img_gray                      (src_gray                       ),
        
        // ���ͼ���� (���ӵ���Ӧ�����)
        .post_img_vsync                    (dst_vsync                      ),
        .post_img_href                     (dst_href                       ),
        .post_img_gray                     (dst_gray                       ) 
    );

    // ��Ӧ�����/ת���� (������Ƶ����д���ļ�)
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
        .i_hsync                           (~dst_href                      ),// ���ӵ�һ�������ֵ, ���ؼ�
        .i_data_valid                      (dst_href                       ),
        .i_data_r                          (dst_gray                       ),
        .i_data_g                          (dst_gray                       ),
        .i_data_b                          (dst_gray                       ),
        .o_frame_count                     (frames_dumped_count            ) 
    );

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // ������������� (Simulation Control & Tasks)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // ������������֡��ת����ֹͣ����
    always @(posedge clk_in2) begin
        if (frames_dumped_count >= FRAMES_TO_SIM) begin
            $display("-------------------------------------------");
            $display("[%0t ns] ����ɹ�������", $time);
            $display("--- ��ת�� %0d ֡ͼ�� ---", FRAMES_TO_SIM);
            $display("--- ����ļ�: %s_%0d.txt ---", OUTPUT_FILE_BASE, FRAMES_TO_SIM-1);
            $display("-------------------------------------------");
            $finish;
        end
    end
    
    // ���һ����ʱ����, ��ֹ����������������Զ������ȥ
    initial begin
        #(CLK1_PERIOD * (SRC_IMG_WIDTH + SRC_H_BLANK) * (SRC_IMG_HEIGHT + SRC_V_BLANK) * (FRAMES_TO_SIM + 1));
        $display("-------------------------------------------");
        $error("[%0t ns] ���泬ʱ����������ʱ�������", $time);
        $display("-------------------------------------------");
        $finish;
    end

    // �������ɲ���ͼ���ļ� (����, һ��ˮƽ����ͼ��) ������
    task automatic generate_input_image(string filename, integer width, integer height);
integer                                  file_handle            ;
integer                                  i,                   j;
integer                                  gray_value             ;
        
        file_handle = $fopen(filename, "w");
        if (file_handle == 0) begin
            $error("TB ����: �޷����ļ� '%s' ����д�롣", filename);
            $finish;
        end

        $display("TB ��Ϣ: �������ɲ���ͼ�� '%s' (%0d x %0d)...", filename, width, height);
        for (i = 0; i < height; i = i + 1) begin
            for (j = 0; j < width; j = j + 1) begin
                // һ����0��255�ļ�ˮƽ����
                gray_value = (j * 255) / (width - 1);
                $fwrite(file_handle, "%x\n", gray_value);
            end
        end
        $fclose(file_handle);
        $display("TB ��Ϣ: ����ͼ��������ϡ�");
    endtask

endmodule