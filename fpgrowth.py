class FPTreeNode():

    def __init__(self, name):
        self.name = name
        self.count = 0
        self.parent = None
        self.children = []
        self.next_homonym = None
        self.tail = None

    def add_child(self, child):
        self.children.append(child)

    def find_child(self, name):
        for child in self.children:
            if child.name == name:
                return child
        return None


class FPGrowth():

    def __init__(self, dataset, min_sup=0.0):
        self.dataset = dataset
        self.min_sup = min_sup
        self.freq_L1 = {}  
        self.freq_itemsets = []  
    def __get_frequency(self, trans_records):
        rect = {}
        for line in trans_records:
            for item in line:
                rect[item] = rect.get(item, 0)+1
        return rect

    def build_fptree(self):
        if self.dataset is None:
            return
        
        self.freq_L1 = self.__get_frequency(self.dataset)
        tmp_list = []
        tmp_list.extend(self.freq_L1.keys())
        tmp_list.sort(key=lambda x: self.freq_L1[x])
        tmp_dict = {}
        i = 1
        for item in tmp_list:
            tmp_dict[item] = i
            i += 1
        
        for trans_record in self.dataset:
            trans_record.sort(key=lambda x: tmp_dict[x], reverse=True)
        self.__fpgrowth(self.dataset, [])
        
        self.freq_itemsets.sort(key=lambda x: len(x[0]), reverse=False)

    def __fpgrowth(self, cpb, post_model):
        freq_dict = self.__get_frequency(cpb)
        headers = {}
        data_num = len(self.dataset)
        for key in freq_dict
            
            if freq_dict.get(key) / data_num >= self.min_sup:
                node = FPTreeNode(key)
                node.count = freq_dict[key]
                headers[key] = node
        tree_root = self.__build_subtree(cpb, headers)
        if len(tree_root.children) == 0:
            return
        for header in headers.values():
            rule = []
            rule.append(header.name)
            rule.extend(post_model)
            
            temp = (rule, header.count / data_num)
            self.freq_itemsets.append(temp)
            
            new_post_pattern = []
            new_post_pattern.append(header.name)
            new_post_pattern.extend(post_model)
    
            new_CPB = []
            next_node = header
            while True:
                next_node = next_node.next_homonym
                if next_node is None:
                    break
                count = next_node.count
                
                path = []
                parent = next_node
                while True:
                    parent = parent.parent
                    #null
                    if parent.name is None:
                        break
                    path.append(parent.name)
                path.reverse()
                #counter
                while count > 0:
                    count -= 1
                    new_CPB.append(path)
            self.__fpgrowth(new_CPB, new_post_pattern)

    def __build_subtree(self, trans_records, headers):
        
        root = FPTreeNode(None)
        for trans_record in trans_records:
            record = []
            record.extend(trans_record)
            subtree_root = root
            tmpRoot = None
            if len(root.children) != 0:
                
                while len(record) != 0:
                    tmpRoot = subtree_root.find_child(record[0])
                    if tmpRoot is None:
                        break
                    tmpRoot.count += 1
                    subtree_root = tmpRoot
                    record.pop(0)
            
            self.__add_nodes(subtree_root, record, headers)
        return root

    def __add_nodes(self, ancestor, record, headers):
        while len(record) > 0:
            item = record.pop(0)
            if item in headers:
                leafnode = FPTreeNode(item)
                leafnode.count = 1
                leafnode.parent = ancestor
                ancestor.add_child(leafnode)
                header = headers[item]
                tail = header.tail
                if tail is None:
                    header.next_homonym = leafnode
                else:
                    tail.next_homonym = leafnode
                header.tail = leafnode
                self.__add_nodes(leafnode, record, headers)

