#include "sim.h"

int main()
{
    // 1.初始化阶段
    q->size = 0; // 硬件堆中现有0个包
    long long sum_packets = 0;
    srand(time(0));
    // manage_metadata.Δ1 = 50000.0, manage_metadata.Δ2 = 100000.0, manage_metadata.Δ3 = 150000.0, manage_metadata.Δ4 = 300000.0; // 初始化Δ
    manage_metadata.Δ1 = 300.0, manage_metadata.Δ2 = 700.0, manage_metadata.Δ3 = 1000.0, manage_metadata.Δ4 = 1500.0; // 初始化Δ

    initialize_memory(); // 初始化SRAM和DRAM
    init_flow_weight();
    initialize_flow_states();
    // 打开输出文件并重定向cout
    ofstream outFile("sim_output_display.txt");
    if (!outFile)
    {
        cerr << "无法打开输出文件！" << endl;
        return 1;
    }
    streambuf *coutBuf = cout.rdbuf(); // 保存原cout缓冲区
    // 2. 用户输入
    int num_cycles;
    cout << "Please enter the number of cycles:";
    cin >> num_cycles;
    cout << "Initial Time is 0.00000" << endl
         << endl;
    int num_cycles_bak = num_cycles;
    int print_ram_cycle = 10;
    while (num_cycles--)
    {
        real_time += 1.0;
        update_flow_states(real_time);
        // 1.生成入包数量/出包数量
        int num_enqueue = gen_num_enqueue_packets();
        int num_dequeue = gen_num_dequeue_packets();
        sum_packets += num_dequeue;
        double middle_rank;
        cout << "Cycle" << num_cycles_bak - num_cycles - 1 << ": " << "Enqueue_Num:" << num_enqueue << ' ' << "Dequeue_Num:" << num_dequeue << endl;
        //  1-1.创建测试数据包
        Packet *test_packets = (Packet *)malloc(sizeof(Packet) * num_enqueue);
        for (int i = 0; i < num_enqueue; i++)
        {
            int flow_id = gen_flow_id();
            if (flow_id == -1)
            {
                flow_id = rand() % (9 - 1) + 1;
            }
            double size = (double)((rand() % (300 - 50 + 1)) + 50) + (double)(1.0 / (double)rand());
            test_packets[i].arrival_time = num_cycles;
            test_packets[i].size = size;
            test_packets[i].flow_id = flow_id;
        }
        Packet *packets[num_enqueue];
        for (int i = 0; i < num_enqueue; i++)
        {
            packets[i] = &test_packets[i];
        }
        // 2.根据所选调度算法类型调用，例如wfq-pre-enqueue(内部调用enqueue_ram)
        if (SCHEDULER_TYPE == WFQ)
        {
            wfq_pre_enqueue(packets, num_enqueue, &middle_rank);
        }

        if (SLIDE_WINDOW_TYPE == NORMAL)
        {
            if (max_rank > manage_metadata.Δ1 + 100.0)
            {
                slide_window_normal(max_rank);
            }
        }
        else if (SLIDE_WINDOW_TYPE == GROWTH_FACTOR)
        {
            slide_window_growth(max_rank, middle_rank);
        }
        else if (SLIDE_WINDOW_TYPE == STATICTICS)
        {
            int empty_bak;
            for (int i = 0; i < SRAM_ROW_NUM; i++)
            {
                for (int j = 0; j < SRAM_COL_NUM; j++)
                {
                    if (SRAM[i][j].rank < 1.0)
                        empty_bak += 1.0;
                }
            }
            if ((empty_bak > 500 && num_cycles < 900) || max_rank > manage_metadata.Δ1 + 300.0)
            {
                slide_window_normal(max_rank);
            }
        }

        // 3.出包处理
        //  3-1.wfq-post-dequeue
        //  3-2.dequeue-ram
        if (SCHEDULER_TYPE == WFQ)
        {
            wfq_post_dequeue(num_dequeue);
        }
        // 4.打印调试信息
        // 5.SRAM和DRAM信息out到output.txt
        if (--print_ram_cycle == 0)
        {
            print_ram_cycle = 10;
            // 打印sram/dram
            cout.rdbuf(outFile.rdbuf());
            cout << "Cycle" << (num_cycles_bak - num_cycles - 1) << endl
                 << endl;
            print_MANAGER();
            print_SRAM();
            print_DRAM();
            cout.rdbuf(coutBuf); // 恢复控制台输出
            outFile.flush();     // 确保数据写入文件
            cout << endl;
        }
    }
    cout << "Hit rata is:" << (double)(hit_num) / (double)sum_packets << endl;
    double num_empty;
    for (int i = 0; i < SRAM_ROW_NUM; i++)
    {
        for (int j = 0; j < SRAM_COL_NUM; j++)
        {
            if (SRAM[i][j].rank < 1.0)
                num_empty += 1.0;
        }
    }
    // printf("%.5f\n", num_empty / 1024.0);
    cout << "SRAM rata is:" << 2.0 * (num_empty / 1024.0) << endl;
    return 0;
}
