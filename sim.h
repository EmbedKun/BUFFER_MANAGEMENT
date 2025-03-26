#ifndef SIM_H
#define SIM_H
// PIFO缓存策略仿真器
#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <limits.h>
#include <fstream>
#include <iomanip>
using namespace std;

/********************  静态参数 ​*​********************/
const int PORT_COUNT = 1; // 交换机端口数量

/******************** 可配置参数 ​********************/
#define SLIDE_WINDOW_TYPE NORMAL // 滑动窗口策略
#define QUEUE_PER_PORT 64        // 每个端口的队列数量
#define MAX_QUEUE_DEPTH 1024     // 单个队列最大缓存包数
#define SIMULATION_TIME 10.0     // 仿真时间（秒）
#define PACKET_RATE 2000.00      // 总包到达速率（包/秒）
#define SERVICE_RATE 1.2e6       // 服务速率（包/秒）
#define SCHEDULER_TYPE WFQ       // 调度策略：SP/RR/DRR/WFQ/VC...
#define MAX_PACKETS 1000
#define N 10

#define PKT_MAX_SIZE 300  // 包最大有效载荷
#define SRAM_ROW_NUM 64   // SRAM行数
#define SRAM_COL_NUM 16   // SRAM列数
#define DRAM_ROW_NUM 1024 // DRAM行数
#define DRAM_COL_NUM 16   // DRAM列数
#define A_ZONE_OFFSET 0   // A区域起始地址相对SRAM的偏移
int B_ZONE_OFFSET = 0;    // B区域起始地址相对DRAM的偏移(不固定)
int C_ZONE_OFFSET = 512;  // C区域起始地址相对DRAM的偏移(不固定)
int D_ZONE_OFFSET = 768;  // D区域起始地址相对DRAM的偏移(不固定)

/******************** 数据结构定义 ​********************/
typedef enum
{
    SP,
    EDF,
    FCFS,
    SRPT,
    SJF,
    LSTF,
    RR,
    WRR,
    DRR,
    DWRR,
    WFQ,
    SCFQ,
    STFQ,
    SFQ
} SchedulerType; // 调度策略枚举

typedef enum
{
    STABLE_FLOW,      // 持续活跃流
    CHURN_FLOW,       // 频繁退出流
    BURSTY_FLOW,      // 突发发包流
    INTERMITTENT_FLOW // 间歇性流
} FlowType;           // 流量类型枚举

typedef struct
{
    FlowType type;
    bool is_active;          // 当前是否活跃
    double last_active_time; // 最后活跃时间
    int burst_counter;       // 突发计数器
    int silent_cycles;       // 静默周期计数器
} FlowState;                 // 流状态结构体

typedef enum
{
    NORMAL,
    GROWTH_FACTOR
} Slide_Window_Type;

typedef struct
{
    int source_port = 1;     // 源端口
    int dest_port = 1;       // 目的端口
    double arrival_time;     // 到达时间
    int size;                // 包大小（字节）
    int flow_id;             // 流标识
    double virt_finish;      // 虚拟完成时间,用于wfq
    double rank;             // PIFO-rank
    char data[PKT_MAX_SIZE]; // 有效载荷
    int row_addr;
    int col_addr;
    int flag;
} Packet;

typedef struct
{
    int size;     // 包大小（字节）
    double rank;  // PIFO-rank
    int flag;     // 指示在A/B/C/D区域
    int flow_id;  // 流标识
    int row_addr; // 包在buffer中的行地址
    int col_addr; // 包在buffer中的列地址
} Metadata;       // 在调度器中的metadata

typedef struct
{
    Packet *buffer[MAX_QUEUE_DEPTH]; // 队列缓存
    int head, tail;                  // 队列指针
    double weight;                   // 队列权重（DRR/WFQ用）
    int deficit_counter;             // 赤字计数器（DRR用）
    double vir_finish_time;          // 虚拟完成时间(FQ类用)
    double vir_start_time;           // 虚拟开始时间(FQ类用)
    long total_bytes;                // 总处理字节
    long dropped_packets;            // 丢包数
    double max_delay;                // 最大延迟
    double total_delay;              // 总延迟
} Queue;                             // 队列/Flow

typedef struct
{
    Queue queues[QUEUE_PER_PORT]; // 每个端口的队列
    int current_queue;            // 当前服务队列（RR用）
} Port;

typedef struct
{
    double history_rank[50][4];
    int history_flow[50][4];
    int history_size[50][4];
    // More...
} History; // 记录历史,从而计算Δ阈值

typedef struct
{
    bool e[N];
    int row[N];
    int col[N];
} DRAM_Addr_Table;

typedef struct
{
    Packet *packets[MAX_PACKETS];
    int size = 0;
} PriorityQueue; // 调度器实现:优先队列（最小堆）

typedef struct
{
    int flag = 0;
    int row_addr;
    int col_addr;
} Cache_idx;

typedef struct
{
    int ring_head_row_idx[4];
    int ring_tail_row_idx[4];
    int ring_head_col_idx[4];
    int ring_tail_col_idx[4];
} ring_idx; // ring-buffer,指示buffer各区域有效数据的首尾位置

typedef struct
{
    double Δ1, Δ2, Δ3, Δ4;
    bool if_forward;
    double last_max_rank = Δ3; // 最大rank区域的下限
    int last_empty_zone = 0;   // 0 for none ,1 for B ,2 for C , 3 for D
    ring_idx ring;
    // SRAM内指针地址
    int sram_row_addr_block_a = 0;
    int sram_col_offset_block_a = 0;
    // DRAM区域B内指针地址
    int dram_row_addr_block_b = B_ZONE_OFFSET;
    int dram_col_offset_block_b = 0;
    // DRAM区域C内指针地址
    int dram_row_addr_block_c = C_ZONE_OFFSET;
    int dram_col_offset_block_c = 0;
    // DRAM区域D内指针地址
    int dram_row_addr_block_d = D_ZONE_OFFSET;
    int dram_col_offset_block_d = 0;
    // SRAM-DRAM-B暂存区
    int sram_dram_b_row_addr = 61;
    int sram_dram_b_col_offset = 0;
    // SRAM-DRAM-C暂存区
    int sram_dram_c_row_addr = 62;
    int sram_dram_c_col_offset = 0;
    // SRAM-DRAM-D暂存区
    int sram_dram_d_row_addr = 63;
    int sram_dram_d_col_offset = 0;
} RAM_MANAGER_METADATA; // 缓存管理

typedef struct
{
    uint32_t **bitmap;      // 二维位图数组
    uint32_t rows;          // 总行数
    uint32_t cols;          // 总列数
    uint32_t words_per_row; // 每行需要的32位字数
} Bitmap;                   // 位图

/******************** 全局变量 ​********************/
Port ports[PORT_COUNT];                                            // 所有交换机端口
RAM_MANAGER_METADATA manage_metadata;                              // 缓存管理器
double real_time = 0.0;                                            // 当前真实时间
double virt_time = 0.0;                                            // 当前虚拟时间,初始为0.FQ类调度算法用
long total_packets = 0;                                            // 总到达包数
long total_dropped = 0;                                            // 总丢包数
long total_processed = 0;                                          // 总处理包数
Packet SRAM[SRAM_ROW_NUM][SRAM_COL_NUM];                           // A
Packet DRAM[DRAM_ROW_NUM][DRAM_COL_NUM];                           // B,C,D
PriorityQueue *q = (PriorityQueue *)malloc(sizeof(PriorityQueue)); // 堆指针
double max_rank = 0.0;                                             // 捕获当前出现过的最大rank
Bitmap bm;                                                         // 位图
DRAM_Addr_Table addr_table[DRAM_ROW_NUM][DRAM_COL_NUM];            // 地址表，用于调度器寻址真包
FlowState flow_states[QUEUE_PER_PORT];                             // 每个流的状态
double historical_max_rank = 0.0;                                  // 全局变量记录历史最大rank（用于Δ4计算）
Slide_Window_Type slide_window_type;
/******************** 工具函数 ​********************/
// 生成随机来包数/出包数(从BMW-Tree出包)
int gen_num_enqueue_packets()
{
    return (rand() % (30 - 2)) + 2;
}
int gen_num_dequeue_packets()
{
    return (rand() % (15 - 2)) + 2;
}

// 辅助函数：获取指定类型的随机活跃流
int get_random_flow_of_type(const vector<int> &active_flows, FlowType type)
{
    vector<int> candidates;
    for (int id : active_flows)
    {
        if (flow_states[id].type == type)
        {
            candidates.push_back(id);
        }
    }
    if (candidates.empty())
        return -1;
    return candidates[rand() % candidates.size()];
}

// 生成随机enqueue包所属流
int gen_flow_id()
{
    vector<int> active_flows;

    // 收集活跃流
    for (int i = 0; i < QUEUE_PER_PORT; i++)
    {
        if (flow_states[i].is_active)
        {
            active_flows.push_back(i);
        }
    }

    if (active_flows.empty())
        return -1; // 无活跃流

    // 不同流量的发包概率权重
    int weights[] = {30, 20, 30, 20}; // 稳定流30%，抖动流20%等

    // 按流量类型加权随机选择
    int selected_type = rand() % 100;
    if (selected_type < 30)
        return get_random_flow_of_type(active_flows, STABLE_FLOW);
    else if (selected_type < 50)
        return get_random_flow_of_type(active_flows, CHURN_FLOW);
    else if (selected_type < 80)
        return get_random_flow_of_type(active_flows, BURSTY_FLOW);
    else
        return get_random_flow_of_type(active_flows, INTERMITTENT_FLOW);
}

// 用于可视化SRAM情况,将结果追加写入文件OUT.txt
void print_SRAM()
{
    for (int i = 0; i < SRAM_ROW_NUM; i++)
    {
        cout << "SRAM[" << i << "]:";
        for (int j = 0; j < SRAM_COL_NUM; j++)
        {
            Packet p = SRAM[i][j];
            cout << fixed << setprecision(2) << p.rank << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

// 用于可视化DRAM情况,将结果追加写入文件OUT.txt
void print_DRAM()
{
    for (int i = 0; i < DRAM_ROW_NUM; i++)
    {
        cout << "DRAM[" << i << "]:";
        for (int j = 0; j < DRAM_COL_NUM; j++)
        {
            Packet p = DRAM[i][j];
            cout << fixed << setprecision(2) << p.rank << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

void print_MANAGER()
{
    cout << fixed << setprecision(2);
    cout << "Δ1:" << manage_metadata.Δ1 << "\t Δ2:" << manage_metadata.Δ2
         << "\t Δ3:" << manage_metadata.Δ3 << "\t Δ4:" << manage_metadata.Δ4 << "\t last_empty_zone:" << manage_metadata.last_empty_zone << endl
         << endl;
}

// 用于打印A/B/C/D区域的置空率
void print_vacancy_rate()
{
}

// 用于echo随机生成的入包/流量
void print_enqueue()
{
}

// 用于echo随机生成的出包/流量
void print_dequeue()
{
}

// 用于插入/寻找DRAM寻址表
// 插入函数：返回插入位置的索引，失败返回-1
void insert_addr(DRAM_Addr_Table *table, int row_val, int col_val)
{
    int last_true = -1;

    // 查找最后一个有效位置
    for (int i = 0; i < N; i++)
    {
        if (table->e[i])
            last_true = i;
    }

    // 存在有效元素时
    if (last_true != -1)
    {
        // 从最后一个有效位置之后insert
        if (last_true != N - 1)
        {
            table->e[last_true + 1] = true;
            table->row[last_true + 1] = row_val;
            table->col[last_true + 1] = col_val;
        }
    }
    // 全表为空时
    else
    {
        if (!table->e[0])
        { // 防御性检查
            table->e[0] = true;
            table->row[0] = row_val;
            table->col[0] = col_val;
        }
    }
}
// 查找函数：返回第一个有效元素的索引，找不到返回-1
int find_first_valid(const DRAM_Addr_Table *table)
{
    for (int i = 0; i < N; i++)
    {
        if (table->e[i])
            return i;
    }
    return -1;
}

// 用于初始化流权重
void init_flow_weight()
{
    for (int i = 0; i < PORT_COUNT; i++)
    {
        double sum_weight;
        for (int j = 0; j < QUEUE_PER_PORT; j++)
        {
            ports[i].queues[j].weight = (rand() % (6 - 1)) + 1;
            sum_weight += ports[i].queues[j].weight;
        }
        for (int j = 0; j < QUEUE_PER_PORT; j++)
        {
            ports[i].queues[j].weight /= sum_weight;
        }
    }
}

// 初始化位图表
Bitmap bitmap_init(uint32_t rows, uint32_t cols)
{
    bm.rows = rows;
    bm.cols = cols;
    bm.words_per_row = (cols + 31) / 32; // 计算每行需要多少个32位字

    // 分配二维数组内存
    bm.bitmap = (uint32_t **)malloc(rows * sizeof(uint32_t *));
    for (uint32_t i = 0; i < rows; i++)
    {
        bm.bitmap[i] = (uint32_t *)calloc(bm.words_per_row, sizeof(uint32_t));
    }
    return bm;
}

// 释放位图资源
void bitmap_free(Bitmap *bm)
{
    for (uint32_t i = 0; i < bm->rows; i++)
    {
        free(bm->bitmap[i]);
    }
    free(bm->bitmap);
}

// 查找第一个空闲位置（原子操作版本）
bool bitmap_find_free(Bitmap *bm, uint32_t *row, uint32_t *col)
{
    for (uint32_t r = 0; r < bm->rows; r++)
    {
        for (uint32_t w = 0; w < bm->words_per_row; w++)
        {
            if (bm->bitmap[r][w] != 0xFFFFFFFF)
            { // 存在空闲位
                uint32_t inverted = ~bm->bitmap[r][w];
                uint32_t bit = __builtin_ffs(inverted) - 1; // 使用GCC内置指令找第一个空闲位
                if (bit < 32)
                {
                    *row = r;
                    *col = w * 32 + bit;
                    if (*col >= bm->cols)
                        continue; // 处理最后一字超出情况

                    // 原子操作标记为已用
                    bm->bitmap[r][w] |= (1U << bit);
                    return true;
                }
            }
        }
    }
    return false;
}

// 硬件堆入队包
void swap(Packet **a, Packet **b)
{
    Packet *temp = *a;
    *a = *b;
    *b = temp;
}
void enqueue(Packet *packet)
{
    if (q->size >= MAX_PACKETS)
        return;
    q->packets[q->size] = packet;
    int i = q->size;
    q->size++;
    while (i > 0 && q->packets[(i - 1) / 2]->rank > q->packets[i]->rank)
    {
        swap(&q->packets[(i - 1) / 2], &q->packets[i]);
        i = (i - 1) / 2;
    }
}

// 硬件堆出队包
void heapify(int i)
{
    int smallest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < q->size && q->packets[left]->rank < q->packets[smallest]->rank)
        smallest = left;
    if (right < q->size && q->packets[right]->rank < q->packets[smallest]->rank)
        smallest = right;
    if (smallest != i)
    {
        swap(&q->packets[i], &q->packets[smallest]);
        heapify(smallest);
    }
}

Packet *dequeue()
{
    if (q->size <= 0)
        return NULL;

    Packet *packet = q->packets[0];
    q->size--;
    q->packets[0] = q->packets[q->size];
    heapify(0);
    return packet;
}

// 清零单个Packet
void init_packet(Packet &p)
{
    memset(&p, 0, sizeof(Packet)); // 所有字段初始化为0
}
// 初始化SRAM和DRAM
void initialize_memory()
{
    // 初始化SRAM
    for (int i = 0; i < SRAM_ROW_NUM; i++)
    {
        for (int j = 0; j < SRAM_COL_NUM; j++)
        {
            init_packet(SRAM[i][j]);
        }
    }
    // 初始化DRAM
    for (int i = 0; i < DRAM_ROW_NUM; i++)
    {
        for (int j = 0; j < DRAM_COL_NUM; j++)
        {
            init_packet(DRAM[i][j]);
        }
    }

    // 初始化地址表
}

/******************** 核心仿真逻辑 ​********************/
// 包到达-ram处理
void enqueue_ram(int rank, Packet *packets)
{
    // 1.根据阈值决定区域
    if (rank <= manage_metadata.Δ1)
    {
        packets->flag = 0;
        if (manage_metadata.sram_row_addr_block_a != manage_metadata.sram_dram_b_row_addr)
        {
            packets->row_addr = manage_metadata.sram_row_addr_block_a, packets->col_addr = manage_metadata.sram_col_offset_block_a;
            SRAM[manage_metadata.sram_row_addr_block_a][manage_metadata.sram_col_offset_block_a++] = (*packets);
            if (manage_metadata.sram_col_offset_block_a == SRAM_COL_NUM)
            {
                manage_metadata.sram_col_offset_block_a = 0;
                manage_metadata.sram_row_addr_block_a++;
            }
            // manage_metadata.ring.ring_tail_row_idx[0] = manage_metadata.sram_row_addr_block_a;
            // manage_metadata.ring.ring_tail_col_idx[0] = manage_metadata.sram_col_offset_block_a;
        }
        else
        {
            cout << "Error:Zone-A is full" << endl;
            manage_metadata.sram_col_offset_block_a = 0; // 还没引入bitmap,先粗略做个FIFO
            manage_metadata.sram_row_addr_block_a = 0;
        }
    }
    else if (rank <= manage_metadata.Δ2)
    {
        if (manage_metadata.last_empty_zone == 0)
        {
            packets->flag = 1;
            packets->row_addr = manage_metadata.dram_row_addr_block_b, packets->col_addr = manage_metadata.dram_col_offset_block_b++;
            SRAM[manage_metadata.sram_dram_b_row_addr][manage_metadata.sram_dram_b_col_offset++] = *packets;
            if (manage_metadata.sram_dram_b_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_b_col_offset = 0;
                manage_metadata.dram_col_offset_block_b = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_b][manage_metadata.dram_col_offset_block_b++] = SRAM[manage_metadata.sram_dram_b_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_b++;
                manage_metadata.dram_col_offset_block_b = 0;
            }
            if (manage_metadata.dram_row_addr_block_b >= C_ZONE_OFFSET)
            {
                cout << "Error:Zone-B is full" << endl;
                manage_metadata.dram_row_addr_block_b = B_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_b = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 1)
        {
            packets->flag = 2;
            packets->row_addr = manage_metadata.dram_row_addr_block_c, packets->col_addr = manage_metadata.dram_col_offset_block_c++;
            SRAM[manage_metadata.sram_dram_c_row_addr][manage_metadata.sram_dram_c_col_offset++] = *packets;
            if (manage_metadata.sram_dram_c_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_c_col_offset = 0;
                manage_metadata.dram_col_offset_block_c = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_c][manage_metadata.dram_col_offset_block_c++] = SRAM[manage_metadata.sram_dram_c_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_c++;
                manage_metadata.dram_col_offset_block_c = 0;
            }
            if (manage_metadata.dram_row_addr_block_c >= D_ZONE_OFFSET)
            {
                cout << "Error:Zone-C is full" << endl;
                manage_metadata.dram_row_addr_block_c = C_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_c = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 2)
        {
            packets->flag = 3;
            packets->row_addr = manage_metadata.dram_row_addr_block_d, packets->col_addr = manage_metadata.dram_col_offset_block_d++;
            SRAM[manage_metadata.sram_dram_d_row_addr][manage_metadata.sram_dram_d_col_offset++] = *packets;
            if (manage_metadata.sram_dram_d_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_d_col_offset = 0;
                manage_metadata.dram_col_offset_block_d = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_d][manage_metadata.dram_col_offset_block_d++] = SRAM[manage_metadata.sram_dram_d_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_d++;
                manage_metadata.dram_col_offset_block_d = 0;
            }
            if (manage_metadata.dram_row_addr_block_d >= DRAM_ROW_NUM)
            {
                cout << "Error:Zone-D is full" << endl;
                manage_metadata.dram_row_addr_block_d = D_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_d = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 3)
        {
            packets->flag = 1;
            packets->row_addr = manage_metadata.dram_row_addr_block_b, packets->col_addr = manage_metadata.dram_col_offset_block_b++;
            SRAM[manage_metadata.sram_dram_b_row_addr][manage_metadata.sram_dram_b_col_offset++] = *packets;
            if (manage_metadata.sram_dram_b_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_b_col_offset = 0;
                manage_metadata.dram_col_offset_block_b = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_b][manage_metadata.dram_col_offset_block_b++] = SRAM[manage_metadata.sram_dram_b_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_b++;
                manage_metadata.dram_col_offset_block_b = 0;
            }
            if (manage_metadata.dram_row_addr_block_b >= C_ZONE_OFFSET)
            {
                cout << "Error:Zone-B is full" << endl;
                manage_metadata.dram_row_addr_block_b = B_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_b = 0;
            }
        }
    }
    else if (rank <= manage_metadata.Δ3)
    {
        if (manage_metadata.last_empty_zone == 0)
        {
            packets->flag = 2;
            packets->row_addr = manage_metadata.dram_row_addr_block_c, packets->col_addr = manage_metadata.dram_col_offset_block_c++;
            SRAM[manage_metadata.sram_dram_c_row_addr][manage_metadata.sram_dram_c_col_offset++] = *packets;
            if (manage_metadata.sram_dram_c_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_c_col_offset = 0;
                manage_metadata.dram_col_offset_block_c = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_c][manage_metadata.dram_col_offset_block_c++] = SRAM[manage_metadata.sram_dram_c_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_c++;
                manage_metadata.dram_col_offset_block_c = 0;
            }
            if (manage_metadata.dram_row_addr_block_c >= D_ZONE_OFFSET)
            {
                cout << "Error:Zone-C is full" << endl;
                manage_metadata.dram_row_addr_block_c = C_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_c = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 1)
        {
            packets->flag = 3;
            packets->row_addr = manage_metadata.dram_row_addr_block_d, packets->col_addr = manage_metadata.dram_col_offset_block_d++;
            SRAM[manage_metadata.sram_dram_d_row_addr][manage_metadata.sram_dram_d_col_offset++] = *packets;
            if (manage_metadata.sram_dram_d_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_d_col_offset = 0;
                manage_metadata.dram_col_offset_block_d = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_d][manage_metadata.dram_col_offset_block_d++] = SRAM[manage_metadata.sram_dram_d_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_d++;
                manage_metadata.dram_col_offset_block_d = 0;
            }
            if (manage_metadata.dram_row_addr_block_d >= DRAM_ROW_NUM)
            {
                cout << "Error:Zone-D is full" << endl;
                manage_metadata.dram_row_addr_block_d = D_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_d = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 2)
        {
            packets->flag = 1;
            packets->row_addr = manage_metadata.dram_row_addr_block_b, packets->col_addr = manage_metadata.dram_col_offset_block_b++;
            SRAM[manage_metadata.sram_dram_b_row_addr][manage_metadata.sram_dram_b_col_offset++] = *packets;
            if (manage_metadata.sram_dram_b_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_b_col_offset = 0;
                manage_metadata.dram_col_offset_block_b = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_b][manage_metadata.dram_col_offset_block_b++] = SRAM[manage_metadata.sram_dram_b_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_b++;
                manage_metadata.dram_col_offset_block_b = 0;
            }
            if (manage_metadata.dram_row_addr_block_b >= C_ZONE_OFFSET)
            {
                cout << "Error:Zone-B is full" << endl;
                manage_metadata.dram_row_addr_block_b = B_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_b = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 3)
        {
            packets->flag = 2;
            packets->row_addr = manage_metadata.dram_row_addr_block_c, packets->col_addr = manage_metadata.dram_col_offset_block_c++;
            SRAM[manage_metadata.sram_dram_c_row_addr][manage_metadata.sram_dram_c_col_offset++] = *packets;
            if (manage_metadata.sram_dram_c_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_c_col_offset = 0;
                manage_metadata.dram_col_offset_block_c = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_c][manage_metadata.dram_col_offset_block_c++] = SRAM[manage_metadata.sram_dram_c_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_c++;
                manage_metadata.dram_col_offset_block_c = 0;
            }
            if (manage_metadata.dram_row_addr_block_c >= D_ZONE_OFFSET)
            {
                cout << "Error:Zone-C is full" << endl;
                manage_metadata.dram_row_addr_block_c = C_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_c = 0;
            }
        }
    }
    else if (rank <= manage_metadata.Δ4)
    {
        if (manage_metadata.last_empty_zone == 0)
        {
            packets->flag = 3;
            packets->row_addr = manage_metadata.dram_row_addr_block_d, packets->col_addr = manage_metadata.dram_col_offset_block_d++;
            SRAM[manage_metadata.sram_dram_d_row_addr][manage_metadata.sram_dram_d_col_offset++] = *packets;
            if (manage_metadata.sram_dram_d_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_d_col_offset = 0;
                manage_metadata.dram_col_offset_block_d = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_d][manage_metadata.dram_col_offset_block_d++] = SRAM[manage_metadata.sram_dram_d_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_d++;
                manage_metadata.dram_col_offset_block_d = 0;
            }
            if (manage_metadata.dram_row_addr_block_d >= DRAM_ROW_NUM)
            {
                cout << "Error:Zone-D is full" << endl;
                manage_metadata.dram_row_addr_block_d = D_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_d = 0;
            }
        }
        else if (manage_metadata.last_empty_zone == 1)
        {
            packets->flag = 1;
            packets->row_addr = manage_metadata.dram_row_addr_block_b, packets->col_addr = manage_metadata.dram_col_offset_block_b++;
            SRAM[manage_metadata.sram_dram_b_row_addr][manage_metadata.sram_dram_b_col_offset++] = *packets;
            if (manage_metadata.sram_dram_b_col_offset == SRAM_COL_NUM)
            {
                manage_metadata.sram_dram_b_col_offset = 0;
                manage_metadata.dram_col_offset_block_b = 0;
                for (int i = 0; i < SRAM_COL_NUM; i++)
                {
                    DRAM[manage_metadata.dram_row_addr_block_b][manage_metadata.dram_col_offset_block_b++] = SRAM[manage_metadata.sram_dram_b_row_addr][i];
                }
                manage_metadata.dram_row_addr_block_b++;
                manage_metadata.dram_col_offset_block_b = 0;
            }
            if (manage_metadata.dram_row_addr_block_b >= C_ZONE_OFFSET)
            {
                cout << "Error:Zone-B is full" << endl;
                manage_metadata.dram_row_addr_block_b = B_ZONE_OFFSET;
                manage_metadata.dram_col_offset_block_b = 0;
            }
        }
    }
    else if (manage_metadata.last_empty_zone == 2)
    {
        packets->flag = 2;
        packets->row_addr = manage_metadata.dram_row_addr_block_c, packets->col_addr = manage_metadata.dram_col_offset_block_c++;
        SRAM[manage_metadata.sram_dram_c_row_addr][manage_metadata.sram_dram_c_col_offset++] = *packets;
        if (manage_metadata.sram_dram_c_col_offset == SRAM_COL_NUM)
        {
            manage_metadata.sram_dram_c_col_offset = 0;
            manage_metadata.dram_col_offset_block_c = 0;
            for (int i = 0; i < SRAM_COL_NUM; i++)
            {
                DRAM[manage_metadata.dram_row_addr_block_c][manage_metadata.dram_col_offset_block_c++] = SRAM[manage_metadata.sram_dram_c_row_addr][i];
            }
            manage_metadata.dram_row_addr_block_c++;
            manage_metadata.dram_col_offset_block_c = 0;
        }
        if (manage_metadata.dram_row_addr_block_c >= D_ZONE_OFFSET)
        {
            cout << "Error:Zone-C is full" << endl;
            manage_metadata.dram_row_addr_block_c = C_ZONE_OFFSET;
            manage_metadata.dram_col_offset_block_c = 0;
        }
    }
    else if (manage_metadata.last_empty_zone == 3)
    {
        packets->flag = 3;
        packets->row_addr = manage_metadata.dram_row_addr_block_d, packets->col_addr = manage_metadata.dram_col_offset_block_d++;
        SRAM[manage_metadata.sram_dram_d_row_addr][manage_metadata.sram_dram_d_col_offset++] = *packets;
        if (manage_metadata.sram_dram_d_col_offset == SRAM_COL_NUM)
        {
            manage_metadata.sram_dram_d_col_offset = 0;
            manage_metadata.dram_col_offset_block_d = 0;
            for (int i = 0; i < SRAM_COL_NUM; i++)
            {
                DRAM[manage_metadata.dram_row_addr_block_d][manage_metadata.dram_col_offset_block_d++] = SRAM[manage_metadata.sram_dram_d_row_addr][i];
            }
            manage_metadata.dram_row_addr_block_d++;
            manage_metadata.dram_col_offset_block_d = 0;
        }
        if (manage_metadata.dram_row_addr_block_d >= DRAM_ROW_NUM)
        {
            cout << "Error:Zone-D is full" << endl;
            manage_metadata.dram_row_addr_block_d = D_ZONE_OFFSET;
            manage_metadata.dram_col_offset_block_d = 0;
        }
    }
}

// 滑动窗口

void slide_window_normal(double current_max_rank)
{
    // 1. 阈值向前滑动
    manage_metadata.last_empty_zone = (manage_metadata.last_empty_zone + 1) % 4;
    if (manage_metadata.last_empty_zone == 0)
        manage_metadata.last_empty_zone = 1;
    manage_metadata.Δ1 = manage_metadata.Δ2; // Δ1 ← 旧Δ2
    manage_metadata.Δ2 = manage_metadata.Δ3; // Δ2 ← 旧Δ3
    manage_metadata.Δ3 = manage_metadata.Δ4; // Δ3 ← 旧Δ4
    manage_metadata.Δ4 += 300.0;
    // 将B区域搬空，搬到SRAM
    if (manage_metadata.last_empty_zone == 1)
    {
        manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
        manage_metadata.sram_col_offset_block_a = 0;
        int row_tmp = B_ZONE_OFFSET;
        int col_tmp = 0;
        while (row_tmp <= manage_metadata.dram_row_addr_block_b)
        {
            for (int i = 1; i <= SRAM_COL_NUM; i++)
            {
                int col_tmp_bak = col_tmp;
                insert_addr(&addr_table[row_tmp][col_tmp], manage_metadata.sram_row_addr_block_a, manage_metadata.sram_col_offset_block_a);
                SRAM[manage_metadata.sram_row_addr_block_a][(manage_metadata.sram_col_offset_block_a++) % SRAM_COL_NUM] = DRAM[row_tmp][(col_tmp++) % DRAM_COL_NUM];
                memset(&DRAM[row_tmp][col_tmp_bak], 0, sizeof(Packet));
            }
            manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
            row_tmp = (row_tmp + 1) % C_ZONE_OFFSET;
        }
        manage_metadata.sram_dram_b_col_offset = 0;
        manage_metadata.dram_row_addr_block_b = B_ZONE_OFFSET;
        manage_metadata.dram_col_offset_block_b = 0;
    } // 将C区域搬空，搬到SRAM
    else if (manage_metadata.last_empty_zone == 2)
    {
        manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
        manage_metadata.sram_col_offset_block_a = 0;
        int row_tmp = C_ZONE_OFFSET;
        int col_tmp = 0;
        while (row_tmp <= manage_metadata.dram_row_addr_block_c)
        {
            for (int i = 1; i <= SRAM_COL_NUM; i++)
            {
                int col_tmp_bak = col_tmp;
                insert_addr(&addr_table[row_tmp][col_tmp], manage_metadata.sram_row_addr_block_a, manage_metadata.sram_col_offset_block_a);
                SRAM[manage_metadata.sram_row_addr_block_a][(manage_metadata.sram_col_offset_block_a++) % SRAM_COL_NUM] = DRAM[row_tmp][(col_tmp++) % DRAM_COL_NUM];
                memset(&DRAM[row_tmp][col_tmp_bak], 0, sizeof(Packet));
            }
            manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
            row_tmp = (row_tmp + 1) % D_ZONE_OFFSET;
        }
        manage_metadata.sram_dram_c_col_offset = 0;
        manage_metadata.dram_row_addr_block_c = C_ZONE_OFFSET;
        manage_metadata.dram_col_offset_block_c = 0;
    }
    else if (manage_metadata.last_empty_zone == 2) // 将D区域搬空，搬到SRAM
    {
        manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
        manage_metadata.sram_col_offset_block_a = 0;
        int row_tmp = D_ZONE_OFFSET;
        int col_tmp = 0;
        while (row_tmp <= manage_metadata.dram_row_addr_block_d)
        {
            for (int i = 1; i <= SRAM_COL_NUM; i++)
            {
                int col_tmp_bak = col_tmp;
                insert_addr(&addr_table[row_tmp][col_tmp], manage_metadata.sram_row_addr_block_a, manage_metadata.sram_col_offset_block_a);
                SRAM[manage_metadata.sram_row_addr_block_a][(manage_metadata.sram_col_offset_block_a++) % SRAM_COL_NUM] = DRAM[row_tmp][(col_tmp++) % DRAM_COL_NUM];
                memset(&DRAM[row_tmp][col_tmp_bak], 0, sizeof(Packet));
            }
            manage_metadata.sram_row_addr_block_a = (manage_metadata.sram_row_addr_block_a + 1) % manage_metadata.sram_dram_b_row_addr;
            row_tmp = (row_tmp + 1) % DRAM_ROW_NUM;
        }
        manage_metadata.sram_dram_d_col_offset = 0;
        manage_metadata.dram_row_addr_block_d = D_ZONE_OFFSET;
        manage_metadata.dram_col_offset_block_d = 0;
    }
}

void slide_window_growth(double current_max_rank)
{
    // 1. 阈值向前滑动
    manage_metadata.last_empty_zone = (manage_metadata.last_empty_zone + 1) % 4;
    if (manage_metadata.last_empty_zone == 0)
        manage_metadata.last_empty_zone = 1;
    manage_metadata.Δ1 = manage_metadata.Δ2; // Δ1 ← 旧Δ2
    manage_metadata.Δ2 = manage_metadata.Δ3; // Δ2 ← 旧Δ3
    manage_metadata.Δ3 = manage_metadata.Δ4; // Δ3 ← 旧Δ4

    // 2. 智能计算新Δ4（核心逻辑）
    const double min_growth = 1.5;  // 最小增长系数
    const double max_growth = 2.0;  // 最大增长系数
    const double hysteresis = 50.0; // 滞后裕量

    // 计算动态增长系数（基于历史趋势）
    double growth_factor = min_growth;
    if (current_max_rank > historical_max_rank)
    {
        // 当检测到新峰值时，扩展幅度与历史增长趋势相关
        double trend = (current_max_rank - historical_max_rank) / historical_max_rank;
        growth_factor = min(max_growth, min_growth + trend * 2.0);
        historical_max_rank = current_max_rank; // 更新历史峰值
    }

    // 应用增长系数并添加滞后裕量
    manage_metadata.Δ4 = current_max_rank * growth_factor + hysteresis;

    // 3. 保护性检查（确保Δ4 > Δ3）
    if (manage_metadata.Δ4 <= manage_metadata.Δ3)
    {
        manage_metadata.Δ4 = manage_metadata.Δ3 * 1.1; // 强制最小增长10%
    }
}

// 初始化流状态
void initialize_flow_states()
{
    for (int i = 0; i < QUEUE_PER_PORT; i++)
    {
        // 分配流量类型（示例配置）
        if (i < 5)
            flow_states[i].type = STABLE_FLOW;
        else if (i < 10)
            flow_states[i].type = CHURN_FLOW;
        else if (i < 15)
            flow_states[i].type = BURSTY_FLOW;
        else
            flow_states[i].type = INTERMITTENT_FLOW;

        flow_states[i].is_active = true;
        flow_states[i].last_active_time = 0;
        flow_states[i].burst_counter = 0;
        flow_states[i].silent_cycles = 0;
    }
}

// 动态更新流状态
void update_flow_states(double current_time)
{
    for (int i = 0; i < QUEUE_PER_PORT; i++)
    {
        switch (flow_states[i].type)
        {
        // 频繁退出流：每5秒有50%概率切换状态
        case CHURN_FLOW:
        {
            if (current_time - flow_states[i].last_active_time > 5.0)
            {
                flow_states[i].is_active = (rand() % 100 < 50);
                flow_states[i].last_active_time = current_time;
            }
            break;
        }

        // 突发流：活跃3秒后静默7秒
        case BURSTY_FLOW:
        {
            if (flow_states[i].is_active)
            {
                if (current_time - flow_states[i].last_active_time > 3.0)
                {
                    flow_states[i].is_active = false;
                    flow_states[i].silent_cycles = 7;
                }
            }
            else
            {
                flow_states[i].silent_cycles--;
                if (flow_states[i].silent_cycles <= 0)
                {
                    flow_states[i].is_active = true;
                    flow_states[i].last_active_time = current_time;
                }
            }
            break;
        }

        // 间歇性流：每10秒激活一次，持续1秒
        case INTERMITTENT_FLOW:
        {
            if (fmod(current_time, 10.0) < 1.0)
            {
                flow_states[i].is_active = true;
            }
            else
            {
                flow_states[i].is_active = false;
            }
            break;
        }

        // 持续活跃流：始终保持活跃
        case STABLE_FLOW:
        default:
            flow_states[i].is_active = true;
            break;
        }
    }
}

/******************** 调度策略实现 ​********************/
// WFQ实现
void wfq_pre_enqueue(Packet *packets[], int num_packets)
{
    for (int i = 0; i < num_packets; i++)
    {
        // 更新流的虚拟时间
        Queue *flow = &ports[0].queues[packets[i]->flow_id];
        packets[i]->virt_finish = ((flow->vir_finish_time > virt_time) && (flow->vir_finish_time < 999900000.0)) ? (flow->vir_finish_time + ((double)packets[i]->size / (PACKET_RATE * flow->weight)))
                                                                                                                 : (virt_time + ((double)packets[i]->size / (PACKET_RATE * flow->weight)));
        flow->vir_finish_time = packets[i]->virt_finish;
        packets[i]->rank = packets[i]->virt_finish;
        double rank = packets[i]->rank;
        max_rank = max_rank > rank ? max_rank : rank;
        // Cache_idx tmp = enqueue_ram(packets[i]->rank, packets[i]);
        // store in SRAM/DRAM,assign addr to pkt
        enqueue_ram(packets[i]->rank, packets[i]);
        // enter BMW-Tree
        enqueue(packets[i]);
        printf("Time %.5f:\t flow_%d_Enque (Length: %.3f) (virt_finish/rank:%f)\t (weight:%f)\n",
               (double)virt_time, packets[i]->flow_id, (double)packets[i]->size, (double)packets[i]->virt_finish, (double)((flow->weight)));
    }
}

void wfq_post_dequeue(int dequeue_num)
{
    while (dequeue_num--)
    {
        Packet *p = dequeue();
        virt_time += (double)(p->size / (PACKET_RATE * (ports[0].queues[p->flow_id].weight)));
        // dequeue_ram();此处实现dequeue_ram函数
        if (ports[0].queues[p->flow_id].total_bytes) // 流非空
        {
            // 由于包是一次性灌输，所以这里不再做enque操作
        }
        if (p->flag == 0)
        {
            memset(&SRAM[p->row_addr][p->col_addr], 0, sizeof(Packet));
            // 发送data
        }
        else if (p->flag != 0)
        {
            int idx = find_first_valid(&addr_table[p->row_addr][p->col_addr]);
            if (idx == -1)
            {
                if (p->flag == 1)
                    cout << "Error:Output Data is dropped by B" << endl;
                else if (p->flag == 2)
                    cout << "Error:Output Data is dropped by C" << endl;
                else if (p->flag == 3)
                    cout << "Error:Output Data is dropped by D" << endl;
            }
            else
            {
                int cur_row = addr_table[p->row_addr][p->col_addr].row[idx];
                int cur_col = addr_table[p->row_addr][p->col_addr].col[idx];
                // 发送data(这里只是清除)
                memset(&SRAM[cur_row][cur_col], 0, sizeof(Packet));
            }
        }
    }
}

#endif