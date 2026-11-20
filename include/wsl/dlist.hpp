#ifndef WSL_DLIST_HPP
#define WSL_DLIST_HPP

#include <compare>
#include <concepts>
#include <cstddef>
#include <iterator>
#include <list>
#include <limits>
#include <ostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace wsl
{
    template <typename T>
    concept streamable = requires(std::ostream &stream, T const &value) {
        stream << value;
    };

    template <typename T, std::signed_integral Offset = std::ptrdiff_t>
    class dlist
    {
    public:
        using offset_type = Offset;
        using position_type = std::make_unsigned_t<offset_type>;
        using value_type = T;
        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using reference = value_type &;
        using const_reference = value_type const &;
        using node_type = std::pair<offset_type, value_type>;
        using type = std::list<node_type>;

        template <bool IsConst>
        class basic_iterator
        {
            using base_iterator = std::conditional_t<IsConst,
                                                     typename type::const_iterator,
                                                     typename type::iterator>;

            base_iterator m_current{};

            explicit basic_iterator(base_iterator current) : m_current(current) {}

            friend class dlist;
            template <bool>
            friend class basic_iterator;

        public:
            using iterator_category = std::bidirectional_iterator_tag;
            using value_type = dlist::value_type;
            using difference_type = std::ptrdiff_t;
            using reference = std::conditional_t<IsConst, value_type const &, value_type &>;
            using pointer = std::conditional_t<IsConst, value_type const *, value_type *>;

            basic_iterator() = default;

            template <bool OtherConst>
                requires IsConst && (!OtherConst)
            basic_iterator(basic_iterator<OtherConst> const &other)
                : m_current(other.m_current)
            {
            }

            reference operator*() const { return m_current->second; }
            pointer operator->() const { return &m_current->second; }

            basic_iterator &operator++()
            {
                ++m_current;
                return *this;
            }

            basic_iterator operator++(int)
            {
                auto result = *this;
                ++*this;
                return result;
            }

            basic_iterator &operator--()
            {
                --m_current;
                return *this;
            }

            basic_iterator operator--(int)
            {
                auto result = *this;
                --*this;
                return result;
            }

            template <bool OtherConst>
            bool operator==(basic_iterator<OtherConst> const &other) const
            {
                return m_current == other.m_current;
            }

            template <bool OtherConst>
            bool operator!=(basic_iterator<OtherConst> const &other) const
            {
                return !(*this == other);
            }

            base_iterator base() const { return m_current; }
        };

        using iterator = basic_iterator<false>;
        using const_iterator = basic_iterator<true>;
        using reverse_iterator = std::reverse_iterator<iterator>;
        using const_reverse_iterator = std::reverse_iterator<const_iterator>;

        struct with_tail_t
        {
            explicit with_tail_t() = default;
        };

        inline static constexpr with_tail_t with_tail{};

    private:
        type m_data{};
        offset_type m_tail{};

        offset_type checked_position(position_type position) const
        {
            if (position > extent() ||
                position > static_cast<position_type>(std::numeric_limits<offset_type>::max()))
                throw std::out_of_range("position exceeds the dlist extent");
            return static_cast<offset_type>(position);
        }

        offset_type checked_length(size_type length) const
        {
            if (length > static_cast<size_type>(std::numeric_limits<offset_type>::max()))
                throw std::out_of_range("length exceeds the offset range");
            return static_cast<offset_type>(length);
        }

        // Premier nœud de position >= position (décalage cumulé courant : un seul parcours).
        iterator node_at_or_after(offset_type position)
        {
            offset_type at = 0;
            for (auto node = begin(); node != end(); ++node)
            {
                at += node.base()->first;
                if (at >= position)
                    return node;
            }
            return end();
        }

        const_iterator node_at_or_after(offset_type position) const
        {
            offset_type at = 0;
            for (auto node = begin(); node != end(); ++node)
            {
                at += node.base()->first;
                if (at >= position)
                    return node;
            }
            return end();
        }

        // Retire un nœud, son écart passant au suivant (ou à la queue), en O(1).
        iterator erase_node(const_iterator position)
        {
            offset_type const gap = position.base()->first;
            iterator next(m_data.erase(position.base()));
            if (next != end())
                next.base()->first += gap;
            else
                m_tail += gap;
            return next;
        }

    public:
        dlist() = default;
        explicit dlist(value_type value)
        {
            m_data.emplace_back(0, std::move(value));
        }
        dlist(with_tail_t, size_type tail)
        {
            if (tail > static_cast<size_type>(std::numeric_limits<offset_type>::max()))
                throw std::out_of_range("tail exceeds the offset range");
            m_tail = static_cast<offset_type>(tail);
        }
        dlist(dlist const &) = default;
        dlist(dlist &&) noexcept = default;
        dlist &operator=(dlist const &) = default;
        dlist &operator=(dlist &&) noexcept = default;

        size_type size() const noexcept { return m_data.size(); }

        bool empty() const noexcept { return m_data.empty(); }

        iterator begin() { return iterator(m_data.begin()); }
        const_iterator begin() const { return const_iterator(m_data.cbegin()); }
        const_iterator cbegin() const { return const_iterator(m_data.cbegin()); }
        iterator end() { return iterator(m_data.end()); }
        const_iterator end() const { return const_iterator(m_data.cend()); }
        const_iterator cend() const { return const_iterator(m_data.cend()); }
        reverse_iterator rbegin() { return reverse_iterator(end()); }
        const_reverse_iterator rbegin() const { return const_reverse_iterator(end()); }
        const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
        reverse_iterator rend() { return reverse_iterator(begin()); }
        const_reverse_iterator rend() const { return const_reverse_iterator(begin()); }
        const_reverse_iterator crend() const { return const_reverse_iterator(cbegin()); }

        reference front() { return m_data.front().second; }
        const_reference front() const { return m_data.front().second; }
        reference back() { return m_data.back().second; }
        const_reference back() const { return m_data.back().second; }

        void clear() noexcept
        {
            m_data.clear();
            m_tail = 0;
        }

    private:
        offset_type accumulated_offset(const_iterator position) const
        {
            offset_type result = 0;
            for (auto node = begin(); node != position; ++node)
                result += node.base()->first;
            return result;
        }

        offset_type node_offset(const_iterator position) const
        {
            return accumulated_offset(position) + position.base()->first;
        }

    public:
        // Position (absolue) de l'élément désigné, en O(distance depuis le début).
        position_type position(const_iterator node) const
        {
            return static_cast<position_type>(node_offset(node));
        }

        position_type extent() const
        {
            auto result = m_data.empty() ? m_tail : node_offset(std::prev(end())) + m_tail;
            if (result < 0 ||
                static_cast<std::make_unsigned_t<offset_type>>(result) >
                    std::numeric_limits<position_type>::max())
                throw std::overflow_error("dlist extent exceeds the position range");
            return static_cast<position_type>(result);
        }

        template <streamable Value = value_type>
        void dump(std::ostream &stream) const
        {
            stream << "dlist{size=" << size() << ", extent=" << extent() << '}';
            if (empty())
                return;

            auto node = cbegin();
            stream << "\n\xE2\x95\x9E\xE2\x95\x90 head  " << node.base()->first;
            stream << "\n\xE2\x94\x9C\xE2\x94\x80 note  ";
            write_value(stream, *node);

            auto previous_position = node_offset(node);
            ++node;
            for (; node != cend(); ++node)
            {
                auto current_position = node_offset(node);
                stream << "\n\xE2\x95\x9F\xE2\x94\x80 gap   "
                       << current_position - previous_position;
                stream << "\n\xE2\x94\x9C\xE2\x94\x80 note  ";
                write_value(stream, *node);
                previous_position = current_position;
            }

            stream << "\n\xE2\x95\x98\xE2\x95\x90 tail  " << m_tail;
        }

    private:
        static void write_value(std::ostream &stream, value_type const &value)
        {
            if constexpr (std::same_as<value_type, std::string>)
                stream << std::quoted(value);
            else
                stream << value;
        }

        void append_node(offset_type offset, value_type value)
        {
            m_data.emplace_back(offset, std::move(value));
        }

        iterator insert_value(position_type position, value_type value)
        {
            auto insertion_position = checked_position(position);
            auto successor = node_at_or_after(insertion_position);
            auto previous_position = successor == begin() ? 0 : node_offset(std::prev(successor));
            auto gap = insertion_position - previous_position;
            auto successor_gap = successor == end() ? m_tail : successor.base()->first;
            auto inserted = m_data.emplace(successor.base(), gap, std::move(value));
            if (successor != end())
                successor.base()->first = successor_gap - gap;
            else
                m_tail = successor_gap - gap;
            return iterator(inserted);
        }

        template <typename... Args>
        iterator emplace_value(position_type position, Args &&...args)
        {
            return insert_value(position, value_type(std::forward<Args>(args)...));
        }

    private:
        iterator insert_impl(position_type position, dlist other)
        {
            auto insertion_position = checked_position(position);
            auto successor = node_at_or_after(insertion_position);
            auto previous_position = successor == begin() ? 0 : node_offset(std::prev(successor));
            auto gap = insertion_position - previous_position;
            auto successor_gap = successor == end() ? m_tail : successor.base()->first;

            auto has_nodes = !other.empty();
            auto inserted = has_nodes ? other.begin() : successor;
            auto first_offset = has_nodes ? inserted.base()->first : 0;
            m_data.splice(successor.base(), other.m_data);
            if (has_nodes)
                inserted.base()->first = gap + first_offset;

            if (successor != end())
                successor.base()->first = other.m_tail + successor_gap - gap;
            else
                m_tail = other.m_tail + successor_gap - gap;
            return inserted;
        }

    public:
        iterator insert(position_type position, value_type value)
        {
            return insert_value(position, std::move(value));
        }

        template <typename... Args>
        iterator emplace(position_type position, Args &&...args)
        {
            return emplace_value(position, std::forward<Args>(args)...);
        }

        iterator insert(position_type position, dlist const &other)
        {
            return insert_impl(position, other);
        }

        iterator insert(position_type position, dlist &&other)
        {
            return insert_impl(position, std::move(other));
        }

        dlist extract(position_type position, size_type length)
        {
            auto first = checked_position(position);
            auto span = checked_length(length);
            if (position > extent() - static_cast<position_type>(span))
                throw std::out_of_range("erase range exceeds the dlist extent");

            dlist result;
            auto successor = node_at_or_after(first + span);
            auto node = node_at_or_after(first);
            auto current_position = node == end() ? offset_type{0} : node_offset(node);
            auto previous_position = first;
            while (node != successor)                                       // positions cumulées au fil du parcours
            {
                result.append_node(current_position - previous_position, std::move(*node));
                previous_position = current_position;
                ++node;
                if (node != successor)
                    current_position += node.base()->first;
            }
            result.m_tail = span - (previous_position - first);

            node = node_at_or_after(first);
            while (node != successor)
                node = erase_node(node);

            if (successor != end())
                successor.base()->first -= span;
            else
                m_tail -= span;
            return result;
        }

        dlist erase(position_type position, size_type length)
        {
            return extract(position, length);
        }

        iterator erase(const_iterator position)
        {
            return erase_node(position);
        }

        iterator erase(const_iterator first, const_iterator last)
        {
            iterator next(m_data.erase(first.base(), first.base()));        // first, rendu modifiable
            while (const_iterator(next) != last)
                next = erase_node(next);
            return next;
        }
    };

    template <typename T, std::signed_integral Offset>
        requires streamable<T>
    std::ostream &operator<<(std::ostream &stream, dlist<T, Offset> const &values)
    {
        values.dump(stream);
        return stream;
    }
}

#endif // WSL_DLIST_HPP