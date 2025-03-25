#include "sim.h"

int main()
{
    // 1.初始化阶段
    q->size = 0; // 硬件堆中现有0个包
    srand(time(0));
    manage_metadata.Δ1 = 50.0, manage_metadata.Δ2 = 100.0, manage_metadata.Δ3 = 150.0, manage_metadata.Δ4 = 300.0; // 初始化Δ
    initialize_memory();                                                                                           // 初始化SRAM和DRAM
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
        cout << "Cycle" << num_cycles_bak - num_cycles - 1 << ": " << "Enqueue_Num:" << num_enqueue << ' ' << "Dequeue_Num:" << num_dequeue << endl;
        //  1-1.创建测试数据包
        Packet *test_packets = (Packet *)malloc(sizeof(Packet) * num_enqueue);
        for (int i = 0; i < num_enqueue; i++)
        {
            int flow_id = gen_flow_id();
            if (flow_id == -1)
                flow_id = rand() % (9 - 1) + 1;
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
            wfq_pre_enqueue(packets, num_enqueue);
        }
        if (max_rank > manage_metadata.Δ4)
        {
            if (SLIDE_WINDOW_TYPE == NORMAL)
            {
                slide_window_normal(max_rank);
            }
            else if (SLIDE_WINDOW_TYPE == GROWTH_FACTOR)
            {
                slide_window_growth(max_rank);
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
    return 0;
}
